"""
File:         combine_and_plot.py
Created:      2020/03/30
Last Changed: 2020/04/04
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
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.local_settings import LocalSettings
from general.df_utilities import save_dataframe, get_basename
from general.utilities import p_value_to_symbol


class CombineAndPlot:
    """
    Main: this class is the main class that calls all other functionality.
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
        self.perm_pvalues_outfile = settings.get_setting("permuted_pvalues_pickle_filename")

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()

        # Combine the pickle files.
        print("Loading data.", flush=True)
        columns, pvalues_data = self.combine_pickles(self.outdir,
                                                     self.pvalues_outfile,
                                                     columns=True)
        _, perm_pvalues = self.combine_pickles(self.outdir,
                                               self.perm_pvalues_outfile)

        # Create a pandas dataframe from the nested list.
        print("Creating p-values dataframe.", flush=True)
        pvalue_df = self.create_df(pvalues_data, columns)

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

        # Visualise null distribution.
        print("Visualizing distributions.", flush=True)
        self.plot_distributions(perm_pvalues, pvalues, self.outdir)

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
                                   self.outdir, stepsize=0.02, max_val=0.05)

        # Print the time.
        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        run_time_hour, run_time_min = divmod(run_time_min, 60)
        print("[manager]\tfinished in  {} hour(s), {} minute(s) and "
              "{} second(s).".format(int(run_time_hour),
                                     int(run_time_min),
                                     int(run_time_sec)), flush=True)

    @staticmethod
    def combine_pickles(indir, filename, columns=False):
        col_list = None
        data = []
        infiles = glob.glob(os.path.join(indir, filename + "*.pkl"))
        start_indices = []
        for fpath in infiles:
            fpath_si = int(fpath.split(filename)[1].split('.')[0])
            start_indices.append(fpath_si)
        start_indices = sorted(start_indices)

        for i, start_index in enumerate(start_indices):
            fpath = os.path.join(indir, filename + str(start_index) + ".pkl")
            with open(fpath, "rb") as f:
                content = pickle.load(f)
                if columns and i == 0:
                    col_list = content[0]
                data.extend(content[1:])
                print("\tLoaded list: {} with length: {}".format(
                    get_basename(fpath), len(content)))
        return col_list, data

    @staticmethod
    def create_df(data, columns):
        df = pd.DataFrame(data, columns=columns)
        df.set_index(df.columns[0], inplace=True)
        df.sort_index(inplace=True)
        df.set_index("-", inplace=True)
        df = df.T

        return df

    @staticmethod
    def plot_distributions(perm_pvalues, pvalues, outdir):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("whitegrid")
        fig, (ax1, ax2) = plt.subplots(1, 2)
        sns.distplot(perm_pvalues,
                     norm_hist=True,
                     color='cornflowerblue',
                     label='null',
                     kde_kws={'clip': (0.0, 1.0)},
                     ax=ax1)
        sns.distplot(pvalues,
                     norm_hist=True,
                     color='firebrick',
                     label='real',
                     kde_kws={'clip': (0.0, 1.0)},
                     ax=ax2)
        for ax, title in zip([ax1, ax2], ["Null", "Alternative"]):
            ax.text(0.5, 1.02,
                    '{} Distribution'.format(title),
                    fontsize=16, weight='bold', ha='center', va='bottom',
                    transform=ax.transAxes)
            ax.set_xlabel('p-value',
                          fontsize=12,
                          fontweight='bold')
            ax.set_ylabel('density',
                          fontsize=12,
                          fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "pvalues_distributions.png"))
        plt.close()

    @staticmethod
    def create_perm_fdr_df(df, pvalues, perm_pvalues, n_perm):
        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)

        len_max_rank = len(str(len(pvalues)))
        len_max_perm_rank = len(str(len(perm_pvalues)))

        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue)
                perm_rank = bisect_left(perm_pvalues, pvalue)
                if (rank > 0) and (perm_rank > 0):
                    fdr_value = (perm_rank / n_perm) / rank
                else:
                    fdr_value = 0
                # print("P-value: {:.2e}\tRank: {:{}d}\tPerm.Rank: {:{}d}\t"
                #       "FDR-value: {:.2e}".format(pvalue,
                #                                    rank, len_max_rank,
                #                                    perm_rank, len_max_perm_rank,
                #                                    fdr_value))
                fdr_df.iloc[row_index, col_index] = fdr_value
        return fdr_df

    @staticmethod
    def create_bh_fdr_df(df, pvalues):
        m = len(pvalues)
        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        fdr_values = []
        len_max_rank = len(str(len(pvalues)))
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue) + 1
                fdr_value = pvalue * (m / rank)
                if fdr_value > 1:
                    fdr_value = 1
                fdr_values.append((rank, row_index, col_index, fdr_value))
                fdr_df.iloc[row_index, col_index] = fdr_value

        fdr_values = sorted(fdr_values, key=lambda x: x[0], reverse=True)
        prev_fdr_value = None
        for (rank, row_index, col_index, fdr_value) in fdr_values:
            # if prev_fdr_value is not None:
                # print("Rank: {}\tRow index: {}\t Col index: {}\t"
                #       "Prev. FDR-Value: {:.2e}\tFDR-Value: {:.2e}".format(rank,
                #                                                           row_index,
                #                                                           col_index,
                #                                                           prev_fdr_value,
                #                                                           fdr_value))
            if prev_fdr_value is not None and fdr_value > prev_fdr_value:
                #print("Changing {} for {}.".format(fdr_value, prev_fdr_value))
                fdr_df.iloc[row_index, col_index] = prev_fdr_value
                prev_fdr_value = prev_fdr_value
            else:
                prev_fdr_value = fdr_value

        return fdr_df

    def create_zscore_df(self, df):
        new_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                pvalue = df.iloc[row_index, col_index]
                zscore = self.get_z_score(pvalue)
                new_df.iloc[row_index, col_index] = zscore

        return new_df

    @staticmethod
    def get_z_score(p_value):
        if p_value > (1.0 - 1e-16):
            p_value = 1e-16
        if p_value < 1e-323:
            p_value = 1e-323

        # stats.norm.isf((1 - 1e-16)) = -8.209536151601387
        # stats.norm.isf(1e-323) = 38.44939448087599
        return stats.norm.isf(p_value)

    def compare_pvalue_scores(self, pvalue_df, perm_fdr, bh_fdr, outdir,
                              stepsize=0.2, max_val=1.0):
        df = pd.DataFrame({"P-value": pvalue_df.melt()["value"],
                           "Permutation FDR": perm_fdr.melt()["value"],
                           "BH FDR": bh_fdr.melt()["value"]})

        # Filter the values of interest.
        indices = []
        for index, row in df.iterrows():
            add = (row <= max_val).any()
            indices.append(add)
        df = df.loc[indices, :]

        # Calculate the lower and upper bound.
        precision = 10 ^ len(str(stepsize).split(".")[1])
        xy_min = math.floor(df.values.min() * precision) / precision - stepsize
        xy_max = math.ceil(df.values.max() * precision) / precision + stepsize

        # Create the plot.
        sns.set()
        sns.set_style("whitegrid")
        sns.set_context("paper", rc={"axes.labelsize": 22,
                                     "axes.fontweight": 'bold',
                                     'xtick.labelsize': 16,
                                     'ytick.labelsize': 16})
        g = sns.pairplot(df, kind='reg', diag_kind='hist', corner=True,
                         height=5, aspect=1,
                         plot_kws={'scatter_kws':
                                       {'facecolors': '#FFFFFF',
                                        'edgecolor': '#FFFFFF'},
                                   'line_kws': {"color": "#FFFFFF"}
                                   },
                         diag_kws={'bins': 25, 'color': "#606060"}
                         )
        # Change the colors.
        g.map_lower(self.my_scatter)
        g.map_lower(self.diagonal_line)
        g.map_lower(self.correlation)
        g.fig.suptitle('Comparing P-values and FDR-values', y=1.05,
                       fontsize=30, fontweight='bold')

        ticks = [round(i / 100, 2) for i in range(int(xy_min*100), int(xy_max*100), int(stepsize*100))]
        for ax in g.axes.flat:
            if ax:
                ax.set_xticks(ticks)
                ax.set_xlim(xy_min, xy_max)
                ax.set_xticklabels(ticks, rotation=45)
                ax.set_yticks(ticks)
                ax.set_ylim(xy_min, xy_max)

        #plt.tight_layout()
        g.savefig(os.path.join(outdir, "pvalue_vs_fdr_comparison_{:.2f}_{:.2f}.png".format(xy_min, xy_max)))

    @staticmethod
    def create_color_map():
        """
        """
        palette = list(Color("#D3D3D3").range_to(Color("#000000"), 1251))
        colors = [x.rgb for x in palette]
        values = [x / 100000 for x in list(range(0, 1251))]
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return value_color_map

    def my_scatter(self, x, y, **kwargs):
        colormap = self.create_color_map()
        diff = round((x - y) ** 2, 5)
        colors = diff.map(colormap).fillna("#000000")
        kwargs['scatter_kws'] = {'facecolors': colors, 'edgecolors': colors}
        kwargs['line_kws'] = {"color": "#D7191C"}
        sns.regplot(x, y, **kwargs)

    @staticmethod
    def diagonal_line(*args, **kwargs):
        plt.plot([-1, 2], [-1, 2], ls="--", c=".3")

    @staticmethod
    def correlation(x, y, **kwargs):
        r, p = stats.spearmanr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.2f} [{}]".format(r, p_value_to_symbol(p)),
                    xy=(.1, .95), xycoords=ax.transAxes, color='#D7191C',
                    fontsize=16, fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > N permutations: {}".format(self.n_permutations))
        print(
            "  > Actual P-values output file: {}*.pkl".format(
                self.pvalues_outfile))
        print("  > Permutated P-values output file: {}*.pkl".format(
            self.perm_pvalues_outfile))
        print("")
