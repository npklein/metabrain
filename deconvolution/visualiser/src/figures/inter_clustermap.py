"""
File:         inter_clustermap.py
Created:      2020/03/16
Last Changed: 2020/04/26
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

    def start(self):
        print("Plotting interaction clustermap")
        self.print_arguments()

        # Plot the z-score matrix.
        print("Plotting the z-score interaction matrix.")
        self.plot(self.inter_df, self.outdir, "zscore", vmin=-8.21, vmax=38.45)

        # Plot the t-value matrix.
        print("Plotting t-value interaction matrix.")
        self.plot(self.inter_tvalue_df, self.outdir, "tvalue ")

    @staticmethod
    def plot(df, outdir, outfile_prefix="", vmin=None, vmax=None):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           vmin=vmin, vmax=vmax,
                           yticklabels=True, xticklabels=False,
                           figsize=(12, (.2 * (len(df.index)))))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "{}_clustermap.png".format(outfile_prefix)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction z-score matrix shape: {}".format(self.inter_df.shape))
        print("  > Interaction t-value matrix shape: {}".format(self.inter_tvalue_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
