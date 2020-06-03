"""
File:         inter_eqtl_effect_celltype.py
Created:      2020/04/20
Last Changed: 2020/06/03
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
import os

# Third party imports.
import pandas as pd
import matplotlib
import itertools
matplotlib.use('Agg')
import upsetplot as up
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLEffectCelltype:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_celltype')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42

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
        celltype_mediated_eqtls = set()
        for (prefix, suffix) in self.cellmap_methods:
            data = {}

            for celltype in self.celltypes:
                method_celltype = celltype
                if method_celltype == "Microglia":
                    method_celltype = "Macrophage"

                df = self.eqtlinter_df.loc[:,
                     prefix + method_celltype + suffix].copy()
                eqtls = set(df.loc[df == 1].index.values)
                celltype_mediated_eqtls.update(eqtls)
                data[method_celltype] = eqtls

            # Plot.
            self.upsetplot(data, prefix.replace("_", "") + suffix,
                           self.outdir, self.extension)

        print("Plotting {}".format(self.marker_genes))
        name = self.marker_genes.replace("_", "")
        data = {}
        for celltype in self.celltypes:
            df = self.eqtlinter_df.loc[:,
                 self.eqtlinter_df.columns.str.startswith(
                     self.marker_genes + celltype + "_")].copy()
            best_marker = df.sum(axis=0).idxmax()
            eqtls = set(df.loc[df[best_marker] == 1].index.values)
            data[celltype] = eqtls

            eqtls = set(df.loc[df.sum(axis=1) > 0].index.values)
            celltype_mediated_eqtls.update(eqtls)

        self.upsetplot(data, name, self.outdir, self.extension)

        print("Plotting all eQTLs.")
        total = len(set(self.eqtlinter_df.index.values))
        part = len(celltype_mediated_eqtls)
        self.plot_pie(total, part, self.outdir, self.extension)

        print("Plotting {} marker genes separately".format(self.marker_genes))
        for celltype in self.celltypes:
            df = self.eqtlinter_df.loc[:,
                 self.eqtlinter_df.columns.str.startswith(
                     self.marker_genes + celltype + "_")].copy()
            data = {}
            for column in df:
                gene_name = column.split("_")[2]
                eqtls = set(df.loc[df[column] == 1, column].index.values)
                data[gene_name] = eqtls
            self.upsetplot(data, "{}{}".format(self.marker_genes, celltype),
                           self.outdir, self.extension)

        print("Plotting all methods combined")
        celltype_data = {}
        for celltype in self.celltypes:
            data = {}

            # Add the marker genes.
            marker_df = self.eqtlinter_df.loc[:,
                        self.eqtlinter_df.columns.str.startswith(
                            self.marker_genes + celltype + "_")].copy()
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
                df = self.eqtlinter_df.loc[:,
                     prefix + method_celltype + suffix].copy()
                eqtls = set(df.loc[df == 1].index.values)
                data[prefix.replace("_", "") + suffix] = eqtls

            # Plot.
            self.upsetplot(data, title, self.outdir, self.extension)

            all_eQTLs_of_celltype = []
            for value in data.values():
                all_eQTLs_of_celltype.extend(value)
            celltype_data[celltype] = set(all_eQTLs_of_celltype)

        print("Plotting all celltypes combined")
        self.upsetplot(celltype_data, "Cell-Type_eQTLs", self.outdir,
                       self.extension)

    def upsetplot(self, data, title, outdir, extension):
        counts = self.count(data)
        up.plot(counts, sort_by='cardinality', show_counts=True)
        plt.suptitle('{}'.format(title.replace("_", " ")), fontsize=18, fontweight='bold')
        plt.savefig(os.path.join(outdir, "{}_upsetplot.{}".format(title, extension)))
        plt.close()

    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        abbreviations = {"Neuron": "neuro", "Oligodendrocyte": "oligo",
                         "EndothelialCell": "endo", "Microglia": "micro",
                         "Macrophage": "macro", "Astrocyte": "astro",
                         "CellMapNNLS": "NNLS", "CellMapPCA_PC1": "PCA",
                         "CellMapNMF_C1": "NMF", "McKenzie": "McKe"}
        abbr_cols = []
        for col in cols:
            if col in abbreviations.keys():
                abbr_cols.append(abbreviations[col])
            else:
                abbr_cols.append(col)

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=abbr_cols))
        s.name = "value"
        return s

    @staticmethod
    def plot_pie(total, part, outdir, extension):
        labels = ['YES', 'NO']
        sizes = [part, total - part]
        explode = (0.1, 0)

        fig, ax = plt.subplots()
        _, _, autotexts = ax.pie(sizes,
                                 explode=explode,
                                 labels=labels,
                                 colors=["#6495ED", "#808080"],
                                 startangle=90,
                                 autopct=lambda pct: "{:.1f}%\n(n={:d})".format(
                                     pct, int(pct / 100. * sum(sizes))),
                                 textprops=dict(fontweight="bold",
                                                fontsize=20))
        for autotext in autotexts:
            autotext.set_color('white')
        ax.axis('equal')

        fig.suptitle("Cell-Type Mediated eQTLs", fontsize=20, fontweight='bold')
        fig.savefig(os.path.join(outdir, "all_cis-eQTLs.{}".format(extension)))
        plt.show()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL and interaction matrix shape: {}".format(
            self.eqtlinter_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
