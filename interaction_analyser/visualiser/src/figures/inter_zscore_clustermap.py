"""
File:         inter_zscore_clustermap.py
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

# Local application imports.
from src.utilities import prepare_output_dir


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
        self.inter_df = dataset.get_inter_df()

    def start(self):
        print("Plotting interaction matrix z-scores as distribution plot.")
        self.print_arguments()
        self.plot(self.inter_df, self.outdir)

    @staticmethod
    def plot(df, outdir):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=False,
                           figsize=(12, (.2 * (len(df.index)))))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=5))
        g.fig.suptitle('Interaction Z-Scores Matrix')
        g.savefig(os.path.join(outdir, "zscores_heatmap.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
