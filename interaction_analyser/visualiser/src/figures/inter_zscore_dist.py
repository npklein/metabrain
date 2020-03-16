"""
File:         inter_zscores_dist.py
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


class InterZscoreDist:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_zscores_dist')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        self.inter_df = dataset.get_inter_df()

    def start(self):
        print("Plotting interaction matrix z-scores as distribution plot.")
        self.print_arguments()

        # calculate the z-score cutoff.
        z_score_cutoff = st.norm.ppf(
            0.05 / (self.inter_df.shape[0] * self.inter_df.shape[1]) / 2)

        self.plot(self.inter_df, z_score_cutoff, self.outdir)

    @staticmethod
    def plot(df, z_score_cutoff, outdir):
        sns.set(style="ticks", color_codes=True)
        df = df.T
        dfm = df.melt(var_name='columns')
        g = sns.FacetGrid(dfm, col='columns', col_wrap=10, sharex=False,
                          sharey=False)
        g.map(sns.distplot, 'value')
        g.map(plt.axvline, x=z_score_cutoff, ls='--', c='red')
        g.map(plt.axvline, x=-1 * z_score_cutoff, ls='--', c='red')
        g.set_titles('{col_name}')
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "cov_zscore_distributions.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
