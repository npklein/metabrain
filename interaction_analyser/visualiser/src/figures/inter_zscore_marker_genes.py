"""
File:         inter_zscore_marker_genes.py
Created:      2020/03/16
Last Changed:
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
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as st

# Local application imports.
from src.utilities import prepare_output_dir


class InterZscoreMarkerGenes:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_zscore_marker_genes')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        self.inter_df = dataset.get_inter_df()

    def start(self):
        print("Plotting interaction matrix.")
        self.print_arguments()

        # Calculate the z-score cutoff.
        z_score_cutoff = st.norm.ppf(
            0.05 / (self.inter_df.shape[0] * self.inter_df.shape[1]))
        gini_cutoff = 0.7

        # Subset on the marker genes.
        marker_df = self.inter_df.loc[
                    (self.inter_df.index.str.contains("neurons")) |
                    (self.inter_df.index.str.contains("oligodendrocytes")) |
                    (self.inter_df.index.str.contains("endothelialcells")) |
                    (self.inter_df.index.str.contains("microglia")) |
                    (self.inter_df.index.str.contains("astrocytes")),
                    :]

        # Create a gini dataframe grouped by celltype.
        gini_df = marker_df.copy()
        gini_df = gini_df.abs()
        gini_df = gini_df.loc[:, gini_df.max(axis=0) > z_score_cutoff]
        celltype_gene = gini_df.index.str.split("_", n=1, expand=True)
        gini_df["celltype"] = [x[0] for x in celltype_gene]
        gini_df = gini_df.groupby(['celltype'], as_index=True).sum()

        # Calculate the gini impurity.
        gini_values = gini_df.div(gini_df.sum(axis=0)).pow(2)
        marker_df = marker_df.T
        marker_df["gini_impurity"] = 1 - gini_values.sum(axis=0)
        marker_df["eqtl_celltype"] = gini_values.idxmax()
        marker_df = marker_df.sort_values(by=['eqtl_celltype', 'gini_impurity'])

        # Subset the marker df on gini impurity.
        marker_df = marker_df.loc[marker_df["gini_impurity"] <= gini_cutoff, :]

        # Plot.
        colormap = self.get_colormap(marker_df['eqtl_celltype'].unique())
        self.plot_clustermap(marker_df, colormap, self.outdir)
        self.plot_bars(marker_df, colormap, z_score_cutoff, gini_cutoff,
                       self.outdir)

    def get_colormap(self, celltypes):
        colors = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e",
                  "#2ecc71"]
        return dict(zip(celltypes, colors))

    @staticmethod
    def plot_clustermap(df, colormap, outdir):
        annotation = df[['gini_impurity', 'eqtl_celltype']].copy()
        data = df.drop(['gini_impurity', 'eqtl_celltype'], axis=1)
        row_colors = list(annotation['eqtl_celltype'].map(colormap))
        column_colors = [colormap[x.split("_")[0]] for x in data.columns]

        sns.set(color_codes=True)
        g = sns.clustermap(data.T, center=0, cmap="RdBu_r",
                           col_cluster=False, row_cluster=False,
                           row_colors=column_colors, col_colors=row_colors,
                           yticklabels=True, xticklabels=False,
                           figsize=(12, 9))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=5))
        g.fig.suptitle('Interaction Z-Scores Marker Genes Gini Impurity')
        g.savefig(os.path.join(outdir, "eqtl_zscores_marker_genes.png"))
        plt.close()

    @staticmethod
    def plot_bars(df, colormap, z_score_cutoff, gini_cutoff, outdir):
        df = df["eqtl_celltype"].copy().value_counts().to_frame()
        df.reset_index(inplace=True)
        df.columns = ["index", "counts"]
        df["hue"] = df["index"].map(colormap)

        sns.set(rc={'figure.figsize': (12, 9)})
        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        g = sns.barplot(x="index", y="counts", data=df, palette=df["hue"])
        g.text(0.5, 1.07,
               'eQTL Counts per Celltype',
               fontsize=16, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.text(0.5, 1.02,
               'zscore >= {:.2f} && gini score <= {:.2f}'.format(z_score_cutoff,
                                                                 gini_cutoff),
               fontsize=10, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('counts',
                     fontsize=12,
                     fontweight='bold')
        g.set_xlabel('celltype',
                     fontsize=12,
                     fontweight='bold')
        ax.tick_params(labelsize=12)
        ax.set_xticks(range(len(df.index)))
        ax.set_xticklabels(df["index"], rotation=45)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "eqtl_marker_genes_barpolot.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
