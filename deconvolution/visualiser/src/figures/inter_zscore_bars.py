"""
File:         inter_zscores_bars.py
Created:      2020/03/16
Last Changed: 2020/04/17
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
import math
import os

# Third party imports.
import scipy.stats as stats
from colour import Color
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


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

        data = self.get_counts(cutoff_value=0)
        self.create_plots(data, "Positive ", "[>0]")

        # Calculate the z-score cutoff.
        z_score_cutoff = stats.norm.ppf(0.05 / (self.inter_df.shape[0] * self.inter_df.shape[1]) / 2)
        self.get_counts(cutoff_value=z_score_cutoff)
        self.create_plots(data, "Significant ", "[>{:.2f}]".format(z_score_cutoff))

    def get_counts(self, cutoff_value=None):
        df = self.inter_df.copy()
        if cutoff_value is not None:
            df[df < cutoff_value] = 0
        sums = df.pow(2).sum(axis=1).to_frame().reset_index()
        sums.columns = ["index", "counts"]
        sums["color"] = self.create_color_map(sums["counts"])
        sums.sort_values(by=['counts'], ascending=False, inplace=True)

        return sums

    def create_plots(self, data, title_prefix="", title_suffix=""):
        majorticks = 10 ** math.floor(math.log10(data["counts"].max()))
        self.plot(data, self.outdir, major_ticks=majorticks, fontsize=5, title_prefix=title_prefix, title_suffix=title_suffix)
        self.plot(data, self.outdir, top=10, major_ticks=majorticks, fontsize=14, title_prefix=title_prefix, title_suffix=title_suffix)
        self.plot(data, self.outdir, bottom=10, major_ticks=majorticks, fontsize=14, title_prefix=title_prefix, title_suffix=title_suffix)

    @staticmethod
    def create_color_map(values, size=101):
        """
        """
        palette = list(Color("#8ABBDB").range_to(Color("#344A5A"), size))
        colors = [str(x).upper() for x in palette]

        start = 0
        stop = math.ceil(max(values))
        step = (stop - start) / size

        color_list = []
        for value in values:
            index = int(value / step) - 1
            color_list.append(colors[index])
        return color_list

    @staticmethod
    def plot(df, outdir, top=None, bottom=None, major_ticks=10, fontsize=10,
             title_prefix="", title_suffix=""):
        # Subset the data.
        title_appendix = ''
        file_appendix = ''
        if top is not None:
            df = df.head(top)
            title_appendix = '(Top {})'.format(top)
            file_appendix = '_top{}'.format(top)
        if bottom is not None:
            df = df.tail(bottom)
            title_appendix = '(Bottom {})'.format(bottom)
            file_appendix = '_bottom{}'.format(bottom)

        # Plot.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        minor_ticks = int(major_ticks / 2)
        for i in range(0, int(max(df["counts"])) + (1 * major_ticks), minor_ticks):
            alpha = 0.025
            if i % major_ticks == 0:
                alpha = 0.15
            ax.axvline(i, ls='-', color="#000000", alpha=alpha, zorder=-1)

        g = sns.barplot(x="counts", y="index", data=df,
                        palette=df["color"], orient="h")

        g.text(0.5, 1.05,
               '{}Covariate Interactions {} {}'.format(title_prefix, title_appendix, title_suffix),
               fontsize=22, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('covariate',
                     fontsize=18,
                     fontweight='bold')
        g.set_xlabel('sum(z-score^2)',
                     fontsize=18,
                     fontweight='bold')
        ax.tick_params(labelsize=fontsize)
        ax.set_yticks(range(len(df.index)))
        ax.set_yticklabels(df["index"])
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}cov_zscores_barplot{}.png".format(title_prefix.replace(" ", "").lower(),
                                                                      file_appendix)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
