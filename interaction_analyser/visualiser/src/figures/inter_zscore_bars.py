"""
File:         inter_zscores_bars.py
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
from src.utilities import prepare_output_dir


class InterZscoreBars:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_zscore_bars')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.inter_df = dataset.get_inter_df()

    def start(self):
        print("Plotting interaction matrix z-scores as barplot.")
        self.print_arguments()

        df = self.inter_df.pow(2)
        sums = df.sum(axis=1).to_frame().reset_index()
        sums.columns = ["index", "counts"]
        sums.sort_values(by=['counts'], ascending=False, inplace=True)

        self.plot(sums, self.outdir)

    @staticmethod
    def plot(df, outdir, top=10):
        sns.set(rc={'figure.figsize': (12, 9)})
        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        g = sns.barplot(x="counts", y="index", data=df, palette="Blues_d",
                        orient="h")
        g.text(0.5, 1.05,
               'Top Covariates',
               fontsize=16, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.text(0.5, 1.02,
               '',
               fontsize=12, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('covariate',
                     fontsize=12,
                     fontweight='bold')
        g.set_xlabel('sum(z-score^2)',
                     fontsize=12,
                     fontweight='bold')
        ax.tick_params(labelsize=5)
        ax.set_yticks(range(len(df.index)))
        ax.set_yticklabels(df["index"])
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "cov_zscores_barplot.png"))
        plt.close()

        subset = df.iloc[:top, :]
        sns.set()
        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        g = sns.barplot(x="counts", y="index", data=subset, palette="Blues_d",
                        orient="h")
        g.set_title('Top Covariates')
        g.set_ylabel('covariate',
                     fontsize=14,
                     fontweight='bold')
        g.set_xlabel('sum(z-score^2)',
                     fontsize=14,
                     fontweight='bold')
        ax.tick_params(labelsize=14)
        ax.set_yticks(range(len(subset.index)))
        ax.set_yticklabels(subset["index"])
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "cov_zscores_barplot_top.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
