"""
File:         inter_clustermap.py
Created:      2020/03/16
Last Changed: 2020/04/29
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
import scipy.stats as stats
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.
from general.utilities import prepare_output_dir


class InterClusterMap:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_clustermap')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.inter_df = dataset.get_inter_df()
        self.inter_tvalue_df = dataset.get_inter_tval_df()
        self.tech_covs = dataset.get_tech_covs()
        self.colormap = {"technical covariate": "cornflowerblue",
                         "covariate of interest": "firebrick"}

    def start(self):
        print("Plotting interaction clustermap")
        self.print_arguments()

        # Plot the z-score matrix.
        print("Plotting the z-score interaction matrix.")
        self.plot(self.inter_df, self.outdir, self.tech_covs, self.colormap,
                  "zscore", vmin=-8.21, vmax=38.45)
        self.plot_colorbar(-8.21, 38.45, self.outdir, "zscore")

        # Plot the t-value matrix.
        print("Plotting t-value interaction matrix.")
        self.plot(self.inter_tvalue_df, self.outdir, self.tech_covs,
                  self.colormap, "tvalue")
        self.plot_colorbar(self.inter_tvalue_df.values.min(),
                           self.inter_tvalue_df.values.max(),
                           self.outdir,
                           "tvalue")

        # Plot the legend.
        self.plot_legend(self.colormap, self.outdir)

    @staticmethod
    def plot(df, outdir, tech_covs, colormap, outfile_prefix="",
             vmin=None, vmax=None):
        # Create the row colors.
        row_colors = []
        for x in df.index:
            if x in tech_covs:
                row_colors.append(colormap["technical covariate"])
            else:
                row_colors.append(colormap["covariate of interest"])

        # Plot.
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           vmin=vmin, vmax=vmax, row_colors=row_colors,
                           yticklabels=True, xticklabels=False,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, (.2 * (len(df.index)))))
        g.cax.set_visible(False)
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.savefig(os.path.join(outdir, "{}_clustermap.png".format(outfile_prefix)))
        plt.close()

    @staticmethod
    def plot_colorbar(vmin, vmax, outdir, name_prefix):
        a = np.array([[vmin, 0, vmax]])
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(3, 9))
        img = sns.heatmap(a, center=0, vmin=vmin, vmax=vmax, cmap="RdBu_r",
                          ax=ax1, cbar_ax=ax2)
        ax1.remove()
        plt.savefig(os.path.join(outdir, "{}_colorbar.png".format(name_prefix)), bbox_inches='tight')

    @staticmethod
    def plot_legend(colormap, outdir):
        sns.set(rc={'figure.figsize': (3, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        handles = []
        for name, color in colormap.items():
            handles.append(mpatches.Patch(color=color, label=name))
        ax.legend(handles=handles)
        ax.set_axis_off()
        fig.savefig(os.path.join(outdir, "legend.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction z-score matrix shape: {}".format(self.inter_df.shape))
        print("  > Interaction t-value matrix shape: {}".format(self.inter_tvalue_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
