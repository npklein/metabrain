"""
File:         inter_eqtl_effect_celltype.py
Created:      2020/04/20
Last Changed: 2020/05/18
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from venn import venn
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLEffectCelltype:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_celltype')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.eqtlinter_df = dataset.get_eqtl_and_interactions_df()
        self.celltypes = dataset.get_celltypes()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()
        self.z_score_cutoff = dataset.get_significance_cutoff()
        self.colormap = dataset.get_colormap()

    def start(self):
        print("Plotting cell type mediated interaction eQTLs")
        self.print_arguments()

        print("Plotting deconvolution methods")
        celltype_mediated_eqtls = []
        for (prefix, suffix) in self.cellmap_methods:
            data = {}

            for celltype in self.celltypes:
                method_celltype = celltype
                if method_celltype == "Microglia":
                    method_celltype = "Macrophage"

                df = self.eqtlinter_df.loc[:, prefix + method_celltype + suffix].copy()
                eqtls = set(df.loc[df == 1].index.values)
                celltype_mediated_eqtls.extend(eqtls)
                data[method_celltype] = eqtls

            # Plot.
            self.plot_data(data, prefix.replace("_", "") + suffix, None,
                           self.outdir, self.colormap, use_pallette=True)

        print("Plotting marker genes")
        data = {}
        for celltype in self.celltypes:
            df = self.eqtlinter_df.loc[:, self.eqtlinter_df.columns.str.startswith(self.marker_genes + celltype + "_")].copy()
            best_marker = df.sum(axis=0).idxmax()
            eqtls = set(df.loc[df[best_marker] == 1].index.values)
            data[celltype] = eqtls

            eqtls = set(df.loc[df.sum(axis=1) > 0].index.values)
            celltype_mediated_eqtls.extend(eqtls)

        self.plot_data(data, self.marker_genes, None,
                       self.outdir, self.colormap, use_pallette=True)

        print("Plotting marker genes separately")
        for celltype in self.celltypes:
            df = self.eqtlinter_df.loc[:, self.eqtlinter_df.columns.str.startswith(self.marker_genes + celltype + "_")].copy()
            data = {}
            for column in df:
                gene_name = column.split("_")[2]
                eqtls = set(df.loc[df[column] == 1, column].index.values)
                data[gene_name] = eqtls
            self.plot_data(data,
                           "{}{}".format(self.marker_genes, celltype),
                           celltype,
                           self.outdir, self.colormap)

        print("Plotting all methods combined")
        celltype_data = {}
        for celltype in self.celltypes:
            data = {}

            # Add the marker genes.
            marker_df = self.eqtlinter_df.loc[:, self.eqtlinter_df.columns.str.startswith(self.marker_genes + celltype + "_")].copy()
            best_marker = marker_df.sum(axis=0).idxmax()
            eqtls = set(marker_df.loc[marker_df[best_marker] == 1].index.values)
            data[self.marker_genes.replace("_", "")] = eqtls

            method_celltype = celltype
            title = celltype
            if method_celltype == "Microglia":
                method_celltype = "Macrophage"
                title = method_celltype + "_" + celltype

            # Add the deconvolution methods.
            for (prefix, suffix) in self.cellmap_methods:
                df = self.eqtlinter_df.loc[:, prefix + method_celltype + suffix].copy()
                eqtls = set(df.loc[df == 1].index.values)
                data[prefix.replace("_", "") + suffix] = eqtls

            # Plot.
            self.plot_data(data, title, celltype, self.outdir, self.colormap)

            all_eQTLs_of_celltype = []
            for value in data.values():
                all_eQTLs_of_celltype.extend(value)
            celltype_data[celltype] = set(all_eQTLs_of_celltype)

        print("Plotting all celltypes combined")
        self.plot_data(celltype_data, "Cell-Type_eQTLs", None, self.outdir,
                       self.colormap, use_pallette=True)

        print("Plotting all eQTLs")
        total = len(set(self.eqtlinter_df.index.values))
        part = len(set(celltype_mediated_eqtls))
        self.plot_pie(total, part, self.outdir)

    @staticmethod
    def plot_data(data, title, celltype, outdir, colormap, use_pallette=False):
        sums = []
        for key, value in data.items():
            sums.append([key, len(value)])
        df = pd.DataFrame(sums, columns=["name", "value"])

        overlap = pd.DataFrame(np.nan, index=data.keys(), columns=data.keys())
        for i, key1 in enumerate(overlap.index):
            for j, key2 in enumerate(overlap.columns):
                value1 = data[key1]
                value2 = data[key2]

                if len(value1) == 0 or len(value2) == 0:
                    continue

                n = len(set(value1).intersection(set(value2)))
                overlap.iloc[i, j] = round(n / len(value1), 2)

        sns.set(color_codes=True)
        g = sns.clustermap(overlap, center=0, cmap="RdBu_r",
                           row_cluster=False, col_cluster=False,
                           yticklabels=True, xticklabels=True, square=True,
                           vmin=0, vmax=1, annot=overlap, fmt='',
                           annot_kws={"size": 16, "color": "#000000"},
                           figsize=(12, 12))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=20, rotation=0))
        plt.setp(
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),
                                         fontsize=20, rotation=45))
        g.fig.suptitle('eQTL Overlap', fontsize=25, fontweight='bold')
        g.savefig(os.path.join(outdir, "{}_overlap.png".format(title)))
        plt.close()

        color = 'cornflowerblue'
        if celltype in colormap:
            color = colormap[celltype]

        sns.set(rc={'figure.figsize': (18, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={"width_ratios": [0.4, 0.6]})
        plt.subplots_adjust(top=0.9, bottom=0.1, wspace=0.01, hspace=0.01)

        sns.despine(fig=fig, ax=ax1)

        major_ticks = 10 ** (math.floor(math.log10(max(df["value"].max(), 100))))
        minor_ticks = int(major_ticks / 2)
        for i in range(0, int(max(df["value"])) + (1 * major_ticks), minor_ticks):
            alpha = 0.025
            if i % major_ticks == 0:
                alpha = 0.15
            ax1.axvline(i, ls='-', color="#000000", alpha=alpha, zorder=-1)

        if use_pallette:
            g = sns.barplot(x="value", y="name", palette=colormap, orient="h",
                            data=df, ax=ax1)
        else:
            g = sns.barplot(x="value", y="name", color=color, orient="h",
                            data=df, ax=ax1)

        g.set_xlabel('n signif. interactions',
                     fontsize=14,
                     fontweight='bold')
        g.set_ylabel('',
                     fontsize=14,
                     fontweight='bold')
        ax1.tick_params(labelsize=14)
        #ax1.set_xticklabels(ax1.get_xmajorticklabels(), rotation=45)

        if use_pallette:
            colors = []
            for key in data.keys():
                colors.append(colormap[key])
            venn(data, fontsize=12,  cmap=ListedColormap(colors), alpha=0.6, legend_loc="upper right", ax=ax2)
        else:
            venn(data, fontsize=12, legend_loc="upper right", ax=ax2)

        ax1.text(1.25, 1.05,
                 '{}'.format(title.replace("_", " ")),
                 fontsize=22, weight='bold', ha='center', va='bottom',
                 transform=ax1.transAxes)

        fig.savefig(os.path.join(outdir,
                                 "{}.png".format(title)))
        plt.close()

    @staticmethod
    def plot_pie(total, part, outdir):
        labels = ['YES', 'NO']
        sizes = [part, total - part]
        explode = (0.1, 0)

        fig, ax = plt.subplots()
        _, _, autotexts = ax.pie(sizes,
                                 explode=explode,
                                 labels=labels,
                                 colors=["#6495ED", "#808080"],
                                 startangle=90,
                                 autopct=lambda pct: "{:.1f}%\n(n={:d})".format(pct, int(pct / 100. * sum(sizes))),
                                 textprops=dict(fontweight="bold",
                                                fontsize=20))
        for autotext in autotexts:
            autotext.set_color('white')
        ax.axis('equal')

        fig.suptitle("Cell-Type Mediated eQTLs", fontsize=20, fontweight='bold')
        fig.savefig(os.path.join(outdir, "all_cis-eQTLs.png"))
        plt.show()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL and interaction matrix shape: {}".format(self.eqtlinter_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
