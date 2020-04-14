"""
File:         celltype_profile_pc_vs_marker_genes.py
Created:      2020/04/07
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
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class CelltypeProfileVSMarkerGenes:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir,
                                   'celltype_profile_pc_vs_marker_genes')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.inter_df = dataset.get_inter_df()
        self.celltypes = dataset.get_celltypes()

        # Create color map.
        self.color_map = self.get_colormap()

    def get_colormap(self):
        colors = ["#9b59b6", "#3498db", "#e74c3c", "#34495e",
                  "#2ecc71"]
        return dict(zip(self.celltypes, colors))

    def start(self):
        print("Plotting cell type PC1 vs Marker genes.")
        self.print_arguments()

        # Get the celltype pc rows.
        ct_pc_rows = []
        for index in self.inter_df.index:
            if index.endswith("_PC1"):
                ct_pc_rows.append(index)
        ct_profile_pc_df = self.inter_df.loc[ct_pc_rows, :]

        # Get the marker gene rows.
        marker_indices = []
        for index in self.inter_df.index:
            if ("_" in index) and (
                    index.split("_")[1] in self.celltypes):
                marker_indices.append(index)
        marker_df = self.inter_df.loc[marker_indices, :]

        # Plot.
        self.plot(ct_profile_pc_df, marker_df, self.celltypes, self.color_map,
                  self.outdir)

    @staticmethod
    def plot(ct_profile_pc_df, marker_df, celltypes, color_map, outdir):
        """
        """
        # Calculate number of rows / columns.
        ncols = len(celltypes)
        nrows = len(ct_profile_pc_df.index)

        sns.set(rc={'figure.figsize': (12 * ncols, 9 * nrows)})
        sns.set_style("ticks")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=ncols,
                                nrows=nrows)

        translate_dict = {"neurons": "NeuronMean_PC1",
                          "oligodendrocytes": "OligodendrocyteMean_PC1",
                          "endothelialcells": "EndothelialMean_PC1",
                          "microglia": "MacrophageMean_PC1",
                          "astrocytes": "AstrocyteMean_PC1"}

        row_index = 0
        col_index = 0
        for celltype in celltypes:
            # Get the pc of this celltype.
            pc_data = ct_profile_pc_df.loc[translate_dict[celltype], :]

            # Get the color.
            color = color_map[celltype]

            # Get the marker genes of this celltype.
            marker_genes = []
            for index in marker_df.index:
                if index.split("_")[1] == celltype:
                    marker_genes.append(index)

            # Loop over each marker gene.
            for marker_gene in marker_genes:
                # Creat the subplot.
                ax = fig.add_subplot(grid[row_index, col_index])
                sns.despine(fig=fig, ax=ax)

                # Get the marker gene data.
                marker_data = marker_df.loc[marker_gene, :]

                # Calculate the correlation.
                coef, p = stats.spearmanr(marker_data, pc_data)

                # Plot.
                g = sns.regplot(x=marker_data,
                                y=pc_data,
                                scatter_kws={'facecolors': '#000000',
                                             'edgecolor': '#000000',
                                             'alpha': 0.5},
                                line_kws={"color": color},
                                ax=ax
                                )

                # Add the text.
                ax.annotate(
                    'r = {:.2f} [{}]'.format(coef, p_value_to_symbol(p)),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color=color,
                    fontsize=12,
                    fontweight='bold')

                ax.text(0.5, 1.05,
                        'Celltype: {} Marker gene: {} '.format(
                            marker_gene.split("_")[1],
                            marker_gene.split("_")[2]),
                        fontsize=16, weight='bold', ha='center', va='bottom',
                        color=color,
                        transform=ax.transAxes)

                ax.set_ylabel(translate_dict[celltype],
                              fontsize=12,
                              fontweight='bold')
                ax.set_xlabel('{} expression'.format(marker_gene.split("_")[1]),
                              fontsize=12,
                              fontweight='bold')

                ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
                ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

                # Increment indices.
                col_index += 1
                if col_index == ncols:
                    col_index = 0
                    row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "celltype_profile_pc_vs_marker_genes.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
