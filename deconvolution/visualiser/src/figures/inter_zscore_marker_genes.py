"""
File:         inter_zscore_marker_genes.py
Created:      2020/03/16
Last Changed: 2020/03/18
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

# Local application imports.
from general.utilities import prepare_output_dir


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
        self.marker_df = dataset.get_marker_df()
        self.expr_df = dataset.get_expr_df()

    def start(self):
        print("Plotting interaction matrix.")
        self.print_arguments()

        # Prepare.
        colormap = self.get_colormap(self.marker_df['eqtl_celltype'].unique())
        annotation = self.marker_df[['gini_impurity', 'eqtl_celltype']].copy()
        data = self.marker_df.drop(['gini_impurity', 'eqtl_celltype', 'SNPName', 'ProbeName', 'HGNCName'], axis=1)

        # Plot.
        self.plot_heatmap(annotation, data, colormap, self.outdir)
        self.plot_bars(annotation, colormap, self.outdir)

        self.marker_df.index = self.marker_df["HGNCName"]
        for celltype in self.marker_df['eqtl_celltype'].unique():
            subset = self.marker_df.loc[self.marker_df['eqtl_celltype'] == celltype,
                                        self.marker_df.columns.str.contains("oligodendrocytes")]
            # subset.drop(['gini_impurity', 'eqtl_celltype', 'SNPName', 'ProbeName',
            #      'HGNCName'], axis=1, inplace=True)
            self.plot_clustermap(subset, celltype, self.outdir)

    def get_colormap(self, celltypes):
        colors = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e",
                  "#2ecc71"]
        return dict(zip(celltypes, colors))

    @staticmethod
    def plot_heatmap(annotation, data, colormap, outdir):
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
                                         fontsize=10))
        g.fig.suptitle('Interaction Z-Scores Marker Genes Gini Impurity')
        g.savefig(os.path.join(outdir, "eqtl_zscores_marker_genes.png"))
        plt.close()

    @staticmethod
    def plot_bars(df, colormap, outdir):
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
               '',
               fontsize=10, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('counts',
                     fontsize=14,
                     fontweight='bold')
        g.set_xlabel('celltype',
                     fontsize=14,
                     fontweight='bold')
        ax.tick_params(labelsize=14)
        ax.set_xticks(range(len(df.index)))
        ax.set_xticklabels(df["index"], rotation=45)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "eqtl_marker_genes_barpolot.png"))
        plt.close()

    @staticmethod
    def plot_clustermap(df, celltype, outdir):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=True,
                           figsize=(12, (.15 * (len(df.index)))))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=5))
        plt.setp(
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),
                                         fontsize=5))
        g.fig.suptitle('Interaction Z-Scores {}'.format(celltype))
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "zscores_{}_clustermap.png".format(celltype)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.marker_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
