"""
File:         combine_and_plot.py
Created:      2020/03/30
Last Changed: 2020/04/23
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
from __future__ import print_function
from pathlib import Path
from bisect import bisect_left
from colour import Color
import pickle
import math
import glob
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.local_settings import LocalSettings
from general.df_utilities import save_dataframe, get_basename
from general.utilities import p_value_to_symbol


class CombineAndPlot:
    """
    Main: this class combines the output of different jobs of the manager into
        one result.
    """
    def __init__(self, settings_file):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.n_permutations = settings.get_setting("n_permutations")

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        if not os.path.exists(self.outdir):
            print("Error, output directory does not exist.")

        # Get the pickle filenames.
        self.pvalues_outfile = settings.get_setting("actual_pvalues_pickle_filename")
        self.tvalues_outfile = settings.get_setting("actual_tvalues_pickle_filename")
        self.perm_pvalues_outfile = settings.get_setting("permuted_pvalues_pickle_filename")

        # Class variable.
        self.max_col_value = 1001
        self.precision = 3

    def start(self):
        """
        Method to start the combiner.
        """
        print("Starting interaction analyser - combine and plot.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()

        # Combine the pickle files.
        print("Loading data.", flush=True)
        pcolumns, pvalues_data = self.combine_pickles(self.outdir,
                                                      self.pvalues_outfile,
                                                      columns=True)
        tcolumns, tvalues_data = self.combine_pickles(self.outdir,
                                                      self.tvalues_outfile,
                                                      columns=True)
        _, perm_pvalues = self.combine_pickles(self.outdir,
                                               self.perm_pvalues_outfile)

        # Create a pandas dataframe from the nested list.
        print("Creating p-values dataframe.", flush=True)
        pvalue_df = self.create_df(pvalues_data, pcolumns)
        save_dataframe(df=pvalue_df,
                       outpath=os.path.join(self.outdir,
                                            "pvalue_table.txt.gz"),
                       header=True, index=True)

        # Create a pandas dataframe from the nested list.
        print("Creating t-values dataframe.", flush=True)
        tvalue_df = self.create_df(tvalues_data, pcolumns)
        save_dataframe(df=tvalue_df,
                       outpath=os.path.join(self.outdir,
                                            "tvalue_table.txt.gz"),
                       header=True, index=True)

        exit()

        # Create a dataframe with z-scores.
        print("Creating Z-score dataframe.", flush=True)
        zscore_df = self.create_zscore_df(pvalue_df)
        save_dataframe(df=zscore_df,
                       outpath=os.path.join(self.outdir,
                                            "interaction_table.txt.gz"),
                       header=True, index=True)

        # Get the pvalues from the dataframe.
        pvalues = pvalue_df.melt()["value"].to_list()

        # Sort the lists.
        print("Sorting p-values.", flush=True)
        perm_pvalues = sorted(perm_pvalues)
        pvalues = sorted(pvalues)

        # Visualise distributions.
        print("Visualizing distributions.", flush=True)
        self.plot_distributions(perm_pvalues, pvalues, self.n_permutations,
                                self.outdir)

        # Create the FDR dataframes.
        print("Creating permutation FDR dataframe.", flush=True)
        perm_fdr_df = self.create_perm_fdr_df(pvalue_df,
                                              pvalues,
                                              perm_pvalues,
                                              self.n_permutations)
        # Write the output file.
        save_dataframe(df=perm_fdr_df,
                       outpath=os.path.join(self.outdir,
                                            "perm_fdr_table.txt.gz"),
                       header=True, index=True)

        print("Creating Benjamini-Hochberg FDR dataframe.", flush=True)
        bh_fdr_df = self.create_bh_fdr_df(pvalue_df, pvalues)
        save_dataframe(df=bh_fdr_df,
                       outpath=os.path.join(self.outdir,
                                            "bh_fdr_table.txt.gz"),
                       header=True, index=True)

        # Compare the two pvalue scores.
        print("Creating score visualisation.", flush=True)
        self.compare_pvalue_scores(pvalue_df, perm_fdr_df, bh_fdr_df,
                                   self.outdir)
        self.compare_pvalue_scores(pvalue_df, perm_fdr_df, bh_fdr_df,
                                   self.outdir, stepsize=0.1, max_val=0.05,
                                   rescale=False)

        # Print the time.
        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        run_time_hour, run_time_min = divmod(run_time_min, 60)
        print("[manager]\tfinished in  {} hour(s), {} minute(s) and "
              "{} second(s).".format(int(run_time_hour),
                                     int(run_time_min),
                                     int(run_time_sec)), flush=True)

    @staticmethod
    def combine_pickles(indir, filename, columns=False):
        """
        Method for combining the pickled lists into one nested list.

        :param indir: string, the input directory containing the pickle files.
        :param filename: string, the prefix name of the input file.
        :param columns: boolean, whether or not each pickle file has a column.
        :return col_list: list, the columns of the content.
        :return data: list, a nested list with the combined content.
        """
        # Declare variables.
        col_list = None
        data = []

        # Combine the found files.
        for i, fpath in enumerate(glob.glob(os.path.join(indir,
                                                         filename + "*.pkl"))):
            with open(fpath, "rb") as f:
                content = pickle.load(f)
                if columns and i == 0:
                    col_list = content[0]
                data.extend(content[1:])
                print("\tLoaded list: {} with length: {}".format(
                    get_basename(fpath), len(content)))
            f.close()

        return col_list, data

    @staticmethod
    def create_df(data, columns):
        """
        Method for creating a pandas dataframe from a nested list.

        :param data: list, a nested list with the data.
        :param columns: list, the column names.
        :return df: DataFrame, the created pandas dataframe.
        """
        df = pd.DataFrame(data, columns=columns)
        df.set_index(df.columns[0], inplace=True)
        df.sort_index(inplace=True)

        reference = set(np.arange(df.index.min(), df.index.max() + 1))
        missing = list(df.index.symmetric_difference(reference))
        duplicated = list(df.index[df.index.duplicated(keep="first")])
        print("\tMissing indices: {}".format(missing))
        print("\tDuplicate indices: {}".format(duplicated))

        df.set_index("-", inplace=True)
        df = df.T

        return df

    def create_zscore_df(self, df):
        """
        Method for converting a dataframe of p-values to z-scores.

        :param df: DataFrame, a dataframe containing p-values.
        :return zscore_df: DataFrame, a dataframe containing z-scores
        """
        count = 1
        total = df.shape[0] * df.shape[1]

        zscore_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                zscore = self.get_z_score(pvalue)
                zscore_df.iloc[row_index, col_index] = zscore
                count += 1

        return zscore_df

    @staticmethod
    def get_z_score(p_value):
        """
        Method for converting a p-value to a z-score.

        :param p-value: float, the p-value to convert.
        :param z-score: float, the corresponding z-score.
        """
        if p_value > (1.0 - 1e-16):
            p_value = (1.0 - 1e-16)
        if p_value < 1e-323:
            p_value = 1e-323

        # The lower and upper limit of stats.norm.sf
        # stats.norm.isf((1 - 1e-16)) = -8.209536151601387
        # stats.norm.isf(1e-323) = 38.44939448087599
        return stats.norm.isf(p_value)

    @staticmethod
    def plot_distributions(perm_pvalues, pvalues, n_perm, outdir):
        """
        Method for visualizing the distribution of the null and alternative
        p-values.

        :param perm_pvalues: list, the sorted null model p-values.
        :param pvalues: list, the sorted alternative model p-values.
        :param n_perm: int, the number of permutations performed.
        :param outdir: string, the output directory for the image.
        """
        # Create bins.
        bins = np.linspace(0, 1, 50)
        perm_pvalues_bins = pd.cut(perm_pvalues, bins=bins).value_counts().to_frame()
        pvalues_bins = pd.cut(pvalues, bins=bins).value_counts().to_frame()
        df = perm_pvalues_bins.merge(pvalues_bins, left_index=True, right_index=True)
        df.columns = ["perm_pvalues", "pvalues"]
        df = df.divide(df.sum())
        df["sum"] = (df["pvalues"] + df["perm_pvalues"]) / 2
        df["index"] = (bins[:-1] + bins[1:]) / 2

        max_val = df.values.max()

        sns.set(rc={'figure.figsize': (18, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=False)
        for i, ax in enumerate([ax1, ax2, ax3]):
            left = True
            if i == 0:
                left = False
            sns.despine(fig=fig, ax=ax, left=left)
            ax.axes.xaxis.set_visible(False)
            if left:
                ax.axes.yaxis.set_visible(False)

            for i in range(0, 101, 5):
                alpha = 0.025
                if i % 10 == 0:
                    alpha = 0.15
                ax.axhline(i / 100, ls='-', color="#000000",
                           alpha=alpha, zorder=-1)

        sns.barplot(x="index", y="perm_pvalues", color='cornflowerblue', data=df, ax=ax1)
        sns.barplot(x="index", y="pvalues", color='firebrick', data=df, ax=ax2)
        sns.barplot(x="index", y="sum", color='firebrick', data=df, ax=ax3)
        sns.barplot(x="index", y="perm_pvalues", color='cornflowerblue', data=df / 2, ax=ax3)

        for ax, title in zip([ax1, ax2, ax3], ["Null", "Alternative", "Combined"]):
            ax.text(0.5, 1.02,
                    '{} Distribution'.format(title),
                    fontsize=16, weight='bold', ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.set_xlabel('p-value',
                          fontsize=12,
                          fontweight='bold')
            ax.set_ylabel('ratio',
                          fontsize=12,
                          fontweight='bold')
            ax.set_ylim(0, 0.35)

        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "pvalues_distributions.png"))
        plt.close()

    @staticmethod
    def create_perm_fdr_df(df, pvalues, perm_pvalues, n_perm):
        """
        Method for creating the permutation False Discovery Rate dataframe.

        FDR = (permutation rank / number of permutations) / actual rank

        :param df: DataFrame, the alternative p-value dataframe.
        :param perm_pvalues: list, the sorted null model p-values.
        :param pvalues: list, the sorted alternative model p-values.
        :param n_perm: int, the number of permutations performed.
        :return fdr_df: DataFrame, the permutation FDR dataframe.
        """
        count = 1
        total = df.shape[0] * df.shape[1]

        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue)
                perm_rank = bisect_left(perm_pvalues, pvalue)
                if (rank > 0) and (perm_rank > 0):
                    fdr_value = (perm_rank / n_perm) / rank
                    if fdr_value > 1:
                        fdr_value = 1
                else:
                    fdr_value = 0
                fdr_df.iloc[row_index, col_index] = fdr_value
                count += 1

        return fdr_df

    @staticmethod
    def create_bh_fdr_df(df, pvalues):
        """
        Method for creating the Benjamini-Hochberg False Discovery Rate
        dataframe.

        FDR = p-value * (# p-values / rank)

        :param df: DataFrame, the alternative p-value dataframe.
        :param pvalues: list, the sorted alternative model p-values.
        :return fdr_df: DataFrame, the Benjamini-Hochberg FDR dataframe.
        """
        count = 1
        total = df.shape[0] * df.shape[1]

        m = len(pvalues)
        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        fdr_values = []
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue) + 1
                fdr_value = pvalue * (m / rank)
                if fdr_value > 1:
                    fdr_value = 1
                fdr_values.append((rank, row_index, col_index, fdr_value))
                fdr_df.iloc[row_index, col_index] = fdr_value
                count += 1

        # Make sure the BH FDR is a monotome function. This goes through
        # the FDR values backwords and make sure that the next FDR is always
        # lower or equal to the previous.
        fdr_values = sorted(fdr_values, key=lambda x: x[0], reverse=True)
        prev_fdr_value = None
        for (rank, row_index, col_index, fdr_value) in fdr_values:
            if prev_fdr_value is not None and fdr_value > prev_fdr_value:
                fdr_df.iloc[row_index, col_index] = prev_fdr_value
                prev_fdr_value = prev_fdr_value
            else:
                prev_fdr_value = fdr_value

        return fdr_df

    def compare_pvalue_scores(self, pvalue_df, perm_fdr, bh_fdr, outdir,
                              stepsize=0.2, max_val=None, rescale=False):
        """
        Method for comparing the different p-value dataframes.

        :param pvalue_df: DataFrame, the alternative p-value dataframe.
        :param perm_fdr: DataFrame, the permutation FDR dataframe.
        :param bh_fdr: DataFrame, the Benjamini-Hochberg FDR dataframe.
        :param outdir: string, the output directory for the image.
        :param stepsize: float, the stepsize for the axis ticks.
        :param max_val: float, the maximum p-value to include.
        :param rescale: boolean, whether or not to convert to -log10 scale.
        """
        # Combine the dataframes into one.
        df = pd.DataFrame({"P-value": pvalue_df.melt()["value"],
                           "Permutation FDR": perm_fdr.melt()["value"],
                           "BH FDR": bh_fdr.melt()["value"]})

        # Rescale the values if need be.
        if rescale:
            df = -np.log10(df + 1)
            df.columns = ["-log10({} + 1)".format(x) for x in df.columns]

        # Filter the rows that contain at least one value below the max.
        if max_val is not None:
            indices = []
            for index, row in df.iterrows():
                add = (row <= max_val).any()
                indices.append(add)
            df = df.loc[indices, :]

        # Calculate the lower and upper bound of the axis.
        precision = 10 ^ len(str(float(stepsize)).split(".")[1])
        xy_min = math.floor(df.values.min() * precision) / precision - stepsize
        xy_max = math.ceil(df.values.max() * precision) / precision + stepsize

        # Calculate the maximum color value.
        max_col_value = -np.inf
        for i in range(3):
            for j in range(3):
                if i > j:
                    new_max = ((df.iloc[:, i] - df.iloc[:, j]) ** 2).max()
                    if new_max > max_col_value:
                        max_col_value = new_max
        self.max_col_value = int(max_col_value * (10 ** self.precision))

        # Create the plot.
        sns.set()
        sns.set_style("whitegrid")
        sns.set_context("paper", rc={"axes.labelsize": 22,
                                     "axes.fontweight": 'bold',
                                     'xtick.labelsize': 16,
                                     'ytick.labelsize': 16})
        g = sns.pairplot(df, kind='scatter', diag_kind='hist', corner=True,
                         height=5, aspect=1,
                         plot_kws={'color': '#FFFFFF'},
                         diag_kws={'bins': 25, 'color': "#606060"}
                         )

        # Change the colors.
        g.map_lower(self.my_scatter)
        g.map_lower(self.diagonal_line)
        g.map_lower(self.correlation)
        g.fig.suptitle('Comparing P-values and FDR-values', y=1.05,
                       fontsize=30, fontweight='bold')

        # Change the axis ticks.
        ticks = [round(i / 100, 2) for i in range(int(xy_min*100), int(xy_max*100), int(stepsize*100))]
        for ax in g.axes.flat:
            if ax:
                ax.set_xticks(ticks)
                ax.set_xlim(xy_min, xy_max)
                ax.set_xticklabels(ticks, rotation=45)
                ax.set_yticks(ticks)
                ax.set_ylim(xy_min, xy_max)

        #plt.tight_layout()
        g.savefig(os.path.join(outdir, "pvalue_vs_fdr_comparison_"
                                       "{:.2f}_{:.2f}.png".format(xy_min,
                                                                  xy_max)))

    @staticmethod
    def create_color_map(length, precision):
        """
        Method for creating a color gradient for the p-values.
        """
        palette = list(Color("#D3D3D3").range_to(Color("#000000"), length + 1))
        colors = [x.rgb for x in palette]
        values = [x / (10 ** precision) for x in list(range(0, length + 1))]
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return value_color_map

    def my_scatter(self, x, y, **kwargs):
        """
        Method for mapping a scatterplot to a FacetGrid. This method includes
        adding the correct color gradient.

        :param x: Series, the x-axis values.
        :param y: Series, the y-xis values.
        """
        colormap = self.create_color_map(length=self.max_col_value,
                                         precision=self.precision)
        diff = round((x - y) ** 2, self.precision)
        colors = diff.map(colormap).fillna("#FFFFFF")
        plt.scatter(x, y, c=colors)

    @staticmethod
    def diagonal_line(*args, **kwargs):
        """
        Method for mapping a diagonal line to a FacetGrid.
        """
        plt.plot([-1, 2], [-1, 2], ls="--", c=".3")

    @staticmethod
    def correlation(x, y, **kwargs):
        """
        Method for mapping a spearman regression annotation to a FacetGrid.

        :param x: Series, the x-axis values.
        :param y: Series, the y-xis values.
        """
        r, p = stats.spearmanr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.2f} [{}]".format(r, p_value_to_symbol(p)),
                    xy=(.1, .95), xycoords=ax.transAxes, color='#D7191C',
                    fontsize=16, fontweight='bold')

    def print_arguments(self):
        """
        Method for printing the variables of the class.
        """
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > N permutations: {}".format(self.n_permutations))
        print(
            "  > Actual P-values output file: {}*.pkl".format(
                self.pvalues_outfile))
        print("  > Permutated P-values output file: {}*.pkl".format(
            self.perm_pvalues_outfile))
        print("")
