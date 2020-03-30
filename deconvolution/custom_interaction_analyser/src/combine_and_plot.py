"""
File:         combine_and_plot.py
Created:      2020/03/30
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
        pvalue_df = self.create_df(pvalues_data, columns)

        # Get the pvalues from the dataframe.
        pvalues = pvalue_df.melt()["value"].to_list()

        # Visualise null distribution.
        self.plot_distributions(pvalues, perm_pvalues, self.outdir)

        # Create a dataframe with adjusted pvalues.
        fdr_df = self.create_fdr_df(pvalue_df,
                                    pvalues,
                                    perm_pvalues,
                                    self.n_permutations)

        # Compare the two pvalue scores.
        self.compare_pvalue_scores(pvalue_df, fdr_df, self.outdir)

        # Create a dataframe with z-scores.
        zscore_df = self.create_zscore_df(pvalue_df)

        # Write the output file.
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
        for i, fpath in enumerate(glob.glob(inpath + "*.pkl")):
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
    def plot_distributions(pvalues, perm_pvalues, outdir):
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
    def create_fdr_df(df, pvalues, perm_pvalues, n_perm):
        pvalues = sorted(pvalues)
        perm_pvalues = sorted(perm_pvalues)
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

    @staticmethod
    def compare_pvalue_scores(pvalue_df, adj_pvalue_df, outdir):
        df = pd.DataFrame({"default": pvalue_df.melt()["value"],
                           "adjusted": adj_pvalue_df.melt()["value"]})
        print(df)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("whitegrid")
        g = sns.jointplot(x="default", y="adjusted", data=df, kind="reg",
                          joint_kws={'scatter_kws':
                                         {'facecolors': '#000000',
                                          'edgecolor': '#000000',
                                          'alpha': 0.2},
                                     'line_kws': {"color": "#D7191C"}
                                     },
                          )
        g.annotate(stats.spearmanr)
        g.ax_marg_x.set_xlim(-0.1, 1.1)
        g.ax_marg_y.set_ylim(-0.1, 1.1)
        g.fig.suptitle('P-value and permutation based FDR regression',
                       fontsize=12,
                       fontweight='bold')
        g.ax_joint.set_ylabel('Permutation based FDR',
                              fontsize=8,
                              fontweight='bold')
        g.ax_joint.set_xlabel('Default P-value',
                              fontsize=8,
                              fontweight='bold')
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "pvalue_comparison.png"))

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > N permutations: {}".format(self.n_permutations))
        print(
            "  > Actual P-values output file: {}*.pkl".format(self.pvalues_outfile))
        print("  > Permutated P-values output file: {}*.pkl".format(
            self.perm_pvalues_outfile))
        print("")
