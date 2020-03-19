"""
File:         inter_eqtl_effect_marker_vs_comp.py
Created:      2020/03/18
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
from colour import Color
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class IntereQTLEffectMarkerVSComp:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_marker_vs_comp')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.inter_df = dataset.get_inter_df()
        self.celltypes = dataset.get_celltypes()

        # Create color map.
        self.color_map = self.get_colormap()

    def start(self):
        print("Plotting interaction eQTL plots for marker genes versus "
              "principle components.")
        self.print_arguments()

        # Determine the columns of the celltypes in the covaraite table.
        marker_indices = []
        comp_indices = []
        for index in self.inter_df.index:
            if ("_" in index) and (
                    index.split("_")[0] in self.celltypes):
                marker_indices.append(index)
            elif index.startswith("Comp"):
                comp_indices.append(index)

        # Loop over the components.
        components_zscores = self.inter_df.loc[comp_indices, :].copy().T
        markers_zscores = self.inter_df.loc[marker_indices, :].copy().T
        self.inter_df = None
        self.cov_df = None

        corr_df = pd.DataFrame(0.0,
                               index=comp_indices, columns=self.celltypes)
        for component in components_zscores.columns:
            print("Plotting component: {}".format(component))
            corr_df = self.plot_scatter_grid(components_zscores,
                                             markers_zscores, component,
                                             self.celltypes, self.color_map,
                                             corr_df, self.outdir)

        print("Average absolute spearman correlation per celltype per comp:")
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(corr_df)

        # Plot a clustermap of the average correlations.
        self.plot_clustermap(corr_df, self.outdir)

    def get_colormap(self):
        colors = ["#9b59b6", "#3498db", "#e74c3c", "#34495e",
                  "#2ecc71"]
        return dict(zip(self.celltypes, colors))

    @staticmethod
    def plot_scatter_grid(components, markers, component, celltypes, color_map,
                          corr_df,outdir):
        """
        """
        # Calculate number of rows / columns.
        ncols = int(len(celltypes))
        nrows = int(len(markers.columns) / len(celltypes))

        sns.set(rc={'figure.figsize': (12 * ncols, 9 * nrows)})
        sns.set_style("whitegrid")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=ncols,
                                nrows=nrows)

        # Get the component from the dataframe.
        comp_zscores = components.loc[:, component].to_frame()

        # Create average cor dict.
        avg_cor_dict = {}

        # Plot the markers.
        row_index = 0
        col_index = 0
        for marker in markers.columns:
            # Split the name.
            celltype = marker.split("_")[0]
            gene = marker.split("_")[1]

            # Get the color.
            color = color_map[celltype]

            # Creat the subplot.
            ax = fig.add_subplot(grid[row_index, col_index])
            sns.despine(fig=fig, ax=ax)

            # Get the marker zscores.
            marker_zscores = markers.loc[:, marker].to_frame()
            df = marker_zscores.merge(comp_zscores, left_index=True,
                                      right_index=True)
            df.columns = ["marker", "comp"]

            # calculate axis limits.
            min = df.values.min() * 1.1
            max = df.values.max() * 1.5

            # Calculate the correlation.
            coef, p = stats.spearmanr(df["marker"], df["comp"])
            if celltype in avg_cor_dict.keys():
                (total_coef, total_p, n) = avg_cor_dict[celltype]
                avg_cor_dict[celltype] = (total_coef + abs(coef),
                                          total_p + p,
                                          n + 1)
            else:
                avg_cor_dict[celltype] = (abs(coef), p, 1)

            # Plot the data.
            sns.regplot(x="marker", y="comp", data=df,
                        scatter_kws={'facecolors': '#000000',
                                     'edgecolor': '#000000',
                                     'alpha': 0.5},
                        line_kws={"color": color},
                        ax=ax)

            # Add the text.
            ax.set(ylim=(min, max))
            ax.set(xlim=(min, max))
            ax.annotate(
                'r = {:.2f}, p = {:.2e} [{}]'.format(coef, p,
                                                     p_value_to_symbol(p)),
                xy=(0.03, 0.94),
                xycoords=ax.transAxes,
                color=color,
                fontsize=12,
                fontweight='bold')

            ax.text(0.5, 1.02,
                    'Celltype: {} Marker gene: {} '.format(celltype, gene),
                    fontsize=16, weight='bold', ha='center', va='bottom',
                    color=color,
                    transform=ax.transAxes)

            ax.set_ylabel('{} z-scores'.format(component),
                          fontsize=12)
            ax.set_xlabel('{} z-scores'.format(marker),
                          fontsize=12,
                          fontweight='bold')

            # Increment indices.
            col_index += 1
            if col_index == ncols:
                col_index = 0
                row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "z_scores_{}_vs_all_markers.png".format(
                                     component)))
        plt.close()

        # Calculate and print the average correlations.
        for key, (total_coef, total_p, n) in avg_cor_dict.items():
            avg_coef = total_coef / n
            avg_p = total_p / n
            corr_df.at[component, key] = avg_coef

        return corr_df

    @staticmethod
    def plot_clustermap(df, outdir, size=0.15):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           col_cluster=False,
                           vmin=0, vmax=1,
                           yticklabels=True, xticklabels=True,
                           figsize=((2 * len(df.columns)),
                                    (0.25 * len(df.index))))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=12))
        plt.setp(
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),
                                         fontsize=14, rotation=45))
        g.fig.suptitle('Average absolute Spearman correlation')
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "sprmn_corr_clustermap.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
