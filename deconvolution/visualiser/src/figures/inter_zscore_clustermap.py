"""
File:         inter_zscore_clustermap.py
Created:      2020/03/16
Last Changed: 2020/04/20
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

# Local application imports.
from general.utilities import prepare_output_dir


class InterZscoreClusterMap:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_zscore_clustermap')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.inter_df = dataset.get_inter_df()

    def start(self):
        print("Plotting interaction matrix z-scores as cluster map.")
        self.print_arguments()

        # Plot the complete matrix.
        print("Plotting the complete matrix.")
        self.plot(self.inter_df, self.outdir)

        # Plot the significant values of the matrix.
        print("Plotting significant values of the matrix.")
        z_score_cutoff = stats.norm.ppf(0.05 / (self.inter_df.shape[0] * self.inter_df.shape[1]) / 2)
        data = self.filter_data(self.inter_df,
                                lower_cutoff=abs(z_score_cutoff),
                                upper_cutoff=np.inf)
        self.plot(data, self.outdir, "Significant ",
                  "\n[{:.2f} < x < {:.2f}]".format(z_score_cutoff,
                                                   abs(z_score_cutoff)))

    @staticmethod
    def filter_data(df, lower_cutoff=0, upper_cutoff=0):
        data = df.copy()
        data[(data > lower_cutoff) & (data < upper_cutoff)] = 0
        return data

    @staticmethod
    def plot(df, outdir, title_prefix="", title_suffix=""):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           vmin=-8.209536151601387, vmax=38.44939448087599,
                           yticklabels=True, xticklabels=False,
                           figsize=(12, (.2 * (len(df.index)))))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.fig.suptitle('{}Interaction Z-Scores Matrix{}'.format(title_prefix,
                                                                title_suffix))
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "{}zscores_clustermap.png".format(title_prefix.replace(" ", "_").lower())))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
