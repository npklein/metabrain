"""
File:         combine_and_plot.py
Created:      2020/03/30
Last Changed: 2020/03/31
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
import pickle
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
from general.df_utilities import save_dataframe
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
        self.pvalues_outfile = os.path.join(self.outdir, settings.get_setting(
            "actual_pvalues_pickle_filename"))
        self.perm_pvalues_outfile = os.path.join(self.outdir,
                                                 settings.get_setting(
                                                     "permuted_pvalues_pickle_filename"))

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()

        # Combine the pickle files.
        columns, pvalues_data = self.combine_pickles(self.pvalues_outfile,
                                                     columns=True)
        _, perm_pvalues = self.combine_pickles(self.perm_pvalues_outfile)

        # Create a pandas dataframe from the nested list.
        print("Creating p-values dataframe.")
        pvalue_df = self.create_df(pvalues_data, columns)

        # Get the pvalues from the dataframe.
        pvalues = pvalue_df.melt()["value"].to_list()

        # Sort the lists.
        print("Sorting p-values.")
        perm_pvalues = sorted(perm_pvalues)
        pvalues = sorted(pvalues)

        # Visualise null distribution.
        print("Visualizing distributions.")
        self.plot_distributions(perm_pvalues, pvalues, self.outdir)

        # Create the FDR dataframes.
        print("Creating permutation FDR dataframe.")
        perm_fdr_df = self.create_perm_fdr_df(pvalue_df,
                                              pvalues,
                                              perm_pvalues,
                                              self.n_permutations)
        print("Creating Benjamini-Hochberg FDR dataframe.")
        bh_fdr_df = self.create_bh_fdr_df(pvalue_df, pvalues)

        # Compare the two pvalue scores.
        print("Creating score visualisation.")
        self.compare_pvalue_scores(pvalue_df, perm_fdr_df, bh_fdr_df,
                                   self.outdir)

        # Create a dataframe with z-scores.
        print("Creating Z-score dataframe.")
        zscore_df = self.create_zscore_df(pvalue_df)

        # Write the output file.
        print("Saving Z-score dataframe.")
        save_dataframe(df=zscore_df,
                       outpath=os.path.join(self.outdir,
                                            "interaction_table.txt.gz"),
                       header=True, index=True)

        # Stop timer.
        process_time = time.time() - start_time

        # Print the time.
        minutes, seconds = divmod(process_time, 60)
        print("Finished in {:.0f} minutes and "
              "{:.0f} seconds".format(minutes, seconds))

    @staticmethod
    def combine_pickles(inpath, columns=False):
        col_list = None
        data = []
        infiles = sorted(glob.glob(inpath + "*.pkl"))
        for i, fpath in enumerate(infiles):
            print("Loading {}.".format(fpath))
            with open(fpath, "rb") as f:
                content = pickle.load(f)
                if columns and i == 0:
                    col_list = content[0]
                data.extend(content[1:])
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

    def create_bh_fdr_df(self, df, pvalues):
        m = len(pvalues)
        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        len_max_rank = len(str(len(pvalues)))
        fdr_values = []
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

    @staticmethod
    def calc_adjusted_pvalue(pvalue, m, rank):
        fdr_value = pvalue * (m / rank)
        if fdr_value > 1:
            fdr_value = 1

        return fdr_value

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
        if p_value >= 1:
            return np.nan
        if p_value <= 0:
            return np.nan

        return stats.norm.ppf((1 - p_value))

    def compare_pvalue_scores(self, pvalue_df, perm_fdr, bh_fdr, outdir):
        df = pd.DataFrame({"P-value": pvalue_df.melt()["value"],
                           "Permutation FDR": perm_fdr.melt()["value"],
                           "BH FDR": bh_fdr.melt()["value"]})

        sns.set()
        sns.set_style("whitegrid")
        sns.set_context("paper", rc={"axes.labelsize": 22,
                                     "axes.fontweight": 'bold',
                                     'xtick.labelsize': 16,
                                     'ytick.labelsize': 16})
        g = sns.pairplot(df, kind='reg', diag_kind='hist', corner=True,
                         height=5, aspect=1,
                         plot_kws={'scatter_kws':
                                       {'facecolors': '#000000',
                                        'edgecolor': '#000000',
                                        'alpha': 0.2},
                                   'line_kws': {"color": "#D7191C"}
                                   },
                         diag_kws={'bins': 25, 'color': "#000000"}
                         )
        g.map_lower(self.correlation)
        g.fig.suptitle('Comparing P-values and FDR-values', y=1.05,
                       fontsize=30, fontweight='bold')

        ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        for ax in g.axes.flat:
            if ax:
                ax.set_xticks(ticks)
                ax.set_xlim(-0.1, 1.1)
                ax.set_yticks(ticks)
                ax.set_ylim(-0.1, 1.1)

        #plt.tight_layout()
        g.savefig(os.path.join(outdir, "pvalue_vs_fdr_comparison.png"))

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