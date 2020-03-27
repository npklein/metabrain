"""
File:         main.py
Created:      2020/03/23
Last Changed: 2020/03/27
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
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
from bisect import bisect_left
import pickle
import random
import queue
import gzip
import time
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import uniform
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.local_settings import LocalSettings
from general.utilities import prepare_output_dir
from general.df_utilities import save_dataframe
from .workers import process_worker


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, n_eqtls, n_samples, cores):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param n_eqtls: int, the number of eqtls in the input files.
        :param n_samples: int, the number of samples in the input files.
        :param cores: int, the number of cores to use.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.geno_inpath = settings.get_setting("genotye_datafile")
        self.expr_inpath = settings.get_setting("expression_datafile")
        self.cov_inpath = settings.get_setting("covariates_datafile")
        self.tech_covs = settings.get_setting("technical_covariates")
        self.max_end_time = int(time.time()) + settings.get_setting(
            "max_runtime_in_min") * 60
        self.n_permutations = settings.get_setting("n_permutations")

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Count the number of eQTLs / samples.
        if n_eqtls is None:
            n_eqtls = self.count_n_eqtls()
        if n_samples is None:
            n_samples = self.count_n_samples()

        # Cap the number of cores to the maximum.
        host_cores = mp.cpu_count()
        if cores > host_cores:
            cores = host_cores

        # Set the correct values for chunk size, and cores.
        if n_eqtls == 1:
            chunk_size = 1
            cores = 1
        elif (n_eqtls > 1) and math.floor(n_eqtls / cores) == 0:
            chunk_size = n_eqtls
            cores = n_eqtls
        else:
            chunk_size = math.ceil(n_eqtls / cores / 2)
        self.n_cores = cores
        self.n_eqtls = n_eqtls
        self.n_samples = n_samples
        self.chunk_size = chunk_size

    def count_n_eqtls(self):
        print("Finding the number of eQTLs")
        geno_nrows = self.get_number_of_rows(self.geno_inpath)
        expr_nrows = self.get_number_of_rows(self.expr_inpath)
        if geno_nrows != expr_nrows:
            print("Number of rows in genotype / expression file do not match.")
            exit()

        return geno_nrows

    def count_n_samples(self):
        print("Finding the number of samples")
        geno_ncols = self.get_number_of_columns(self.geno_inpath)
        expr_ncols = self.get_number_of_columns(self.expr_inpath)
        if geno_ncols != expr_ncols:
            print("Number of columns in genotype / expression file do "
                  "not match.")
            exit()

        return geno_ncols

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()
        #
        # # Create queues.
        # input_manager = mp.Manager()
        # input_queue = input_manager.Queue()
        # output_manager = mp.Manager()
        # output_queue = output_manager.Queue()
        #
        # # Fill input queue.
        # print("Loading input queue.")
        # self.fill_input_queue(input_queue)
        #
        # # Fire off the workers.
        # processes = []
        # for worker_id in range(self.n_cores):
        #     processes.append(mp.Process(target=process_worker,
        #                                 args=(worker_id,
        #                                       self.cov_inpath,
        #                                       self.geno_inpath,
        #                                       self.expr_inpath,
        #                                       self.tech_covs,
        #                                       input_queue,
        #                                       output_queue,
        #                                       self.max_end_time,
        #                                       self.chunk_size,
        #                                       self.n_permutations)))
        # for proc in processes:
        #     proc.start()
        #
        # # Empty the output queue.
        # columns = None
        # pvalue_data = []
        # perm_pvalues = []
        # eqtls_done = 0
        # total_zscores = self.n_eqtls * (self.n_permutations + 1)
        # last_print_time = time.time() - 60
        # print("Waiting for output from workers.")
        # while eqtls_done != total_zscores:
        #     now_time = time.time()
        #     if self.max_end_time <= now_time:
        #         print("ERROR: max progress time reached")
        #         exit()
        #     if now_time - last_print_time >= 60:
        #         last_print_time = now_time
        #         print("\tReceived data for {}/{} eQTLs "
        #               "[{:.2f}%]".format(eqtls_done,
        #                                  total_zscores,
        #                                  (100 / total_zscores) * eqtls_done))
        #
        #     try:
        #         if output_queue.empty():
        #             time.sleep(0.1)
        #             continue
        #
        #         # Get the new output
        #         (output_type, output) = output_queue.get(True, timeout=1)
        #         new_index = output[0]
        #
        #         # Check the length.
        #         if columns is not None and len(output) != len(columns):
        #             print("Length of new ouput on index '{}' does not "
        #                   "correspond with the remainder of the "
        #                   "data.".format(new_index))
        #             exit()
        #
        #         # Safe the columns.
        #         if new_index == -1:
        #             if columns is None:
        #                 columns = output
        #             else:
        #                 if output != columns:
        #                     print("Columns are not identical over "
        #                           "subprocesses.")
        #         else:
        #             # Safe the output to the right result.
        #             eqtls_done += 1
        #             if output_type == -1:
        #                 pvalue_data.append(output)
        #             else:
        #                 perm_pvalues.extend(output[2:])
        #
        #     except queue.Empty:
        #         time.sleep(1)
        #         continue
        # print("\tReceived data for {}/{} eQTLs "
        #       "[{:.2f}%]".format(eqtls_done,
        #                          total_zscores,
        #                          (100 / total_zscores) * eqtls_done))
        #
        # # Join the processes.
        # for proc in processes:
        #     proc.join()
        #
        # # Create a pandas dataframe from the nested list.
        # pvalue_df = self.create_df(pvalue_data, columns)
        #
        # # Pickle the files.
        # pvalue_df.to_pickle(os.path.join(self.outdir, "pvalue_df.pkl"))
        # with open(os.path.join(self.outdir, "perm_pvalues.pkl"), "wb") as f:
        #     pickle.dump(perm_pvalues, f)
        # f.close()

        # Load the pickled files.
        pvalue_df = pd.read_pickle(os.path.join(self.outdir, "pvalue_df.pkl"))
        with open(os.path.join(self.outdir, "perm_pvalues.pkl"), "rb") as f:
            perm_pvalues = pickle.load(f)
        f.close()

        # Get the pvalues from the dataframe.
        pvalues = sorted(pvalue_df.melt()["value"].to_list())

        # Visualise null distribution.
        self.plot_distributions(pvalues, perm_pvalues, self.outdir)

        # Create a dataframe with adjusted pvalues.
        adj_pvalue_df = self.create_adj_pvalue_df(pvalue_df,
                                                  pvalues,
                                                  perm_pvalues,
                                                  self.n_permutations)

        # Compare the two pvalue scores.
        self.compare_pvalue_scores(pvalue_df, adj_pvalue_df, self.outdir)

        # Create a dataframe with z-scores.
        zscore_df = self.create_zscore_df(adj_pvalue_df)

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
    def get_number_of_rows(inpath):
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                pass
        f.close()
        return i

    @staticmethod
    def get_number_of_columns(inpath):
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                ncols = len(l.decode().split("\t")) - 1
                break
        f.close()
        return ncols

    def fill_input_queue(self, input_queue):
        # Create a list of start indices.
        start_indices = [i for i in range(1, self.n_eqtls + 1, self.chunk_size)]
        if self.n_eqtls == 1 and self.chunk_size == 1:
            start_indices = [1]

        # Create x random shuffles of the column indices.
        random_shuffles = [self.create_random_shuffle() for _ in
                           range(self.n_permutations)]

        # Start filling the input queue.
        njobs = 0
        for start_index in start_indices:
            # print("input_queue.put(({}, -1, {}))".format(start_index, self.get_column_indices()[:5]))
            input_queue.put((start_index, -1, self.get_column_indices()))
            njobs += 1
            for i, random_shuffle in enumerate(random_shuffles):
                # print("input_queue.put(({}, {}, {}))".format(start_index, i, random_shuffle[:5]))
                input_queue.put((start_index, i, random_shuffle))
                njobs += 1

        print("Number of jobs to complete: {}".format(njobs))

    def get_column_indices(self):
        return [x for x in range(self.n_samples)]

    def create_random_shuffle(self):
        column_indices = self.get_column_indices()
        random.shuffle(column_indices)
        return column_indices

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
                     fit=uniform,
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
        fig.savefig(os.path.join(outdir, "perm_pvalues_distributions.png"))
        plt.close()

    @staticmethod
    def create_adj_pvalue_df(df, pvalues, perm_pvalues, n_perm):
        pvalues = sorted(pvalues)
        perm_pvalues = sorted(perm_pvalues)
        new_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in df.index:
            for col_index in df.columns:
                pvalue = df.at[row_index, col_index]
                rank = bisect_left(pvalues, pvalue) + 1
                perm_rank = bisect_left(perm_pvalues, pvalue) + 1
                adj_pvalue = (perm_rank / n_perm) / rank
                print("P-value: {:.4f}\tRank: {}\tPerm.Rank: {}\t"
                      "Adj.P-value: {:.4f}".format(pvalue, rank,
                                               perm_rank, adj_pvalue))
                new_df.at[row_index, col_index] = adj_pvalue

        return new_df

    def create_zscore_df(self, df):
        new_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in df.index:
            for col_index in df.columns:
                pvalue = df.at[row_index, col_index]
                zscore = self.get_z_score(pvalue)
                new_df.at[row_index, col_index] = zscore

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

        coef, p_value = stats.spearmanr(df["default"], df["adjusted"])

        fig, ax = plt.subplots()
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        g = sns.regplot(x="default", y="adjusted", data=df,
                        scatter_kws={'facecolors': '#000000',
                                     'edgecolor': '#000000',
                                     'alpha': 0.5},
                        line_kws={"color": "#D7191C"},
                        )
        g.set_title('Z-score regression (Spearman r = {:.2f}, '
                    'p = {:.2e})'.format(coef, p_value))
        g.set_ylabel('Adjusted pvalue',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('Default pvalue',
                     fontsize=8,
                     fontweight='bold')
        g.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        g.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        fig.savefig(os.path.join(outdir, "pvalue_comparison.png"))

    def print_arguments(self):
        end_time_string = datetime.fromtimestamp(self.max_end_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        print("Arguments:")
        print("  > Genotype datafile: {}".format(self.geno_inpath))
        print("  > Expression datafile: {}".format(self.expr_inpath))
        print("  > Covariates datafile: {}".format(self.cov_inpath))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > N eQTLs: {}".format(self.n_eqtls))
        print("  > N samples: {}".format(self.n_samples))
        print("  > N cores: {}".format(self.n_cores))
        print("  > Chunk size: {}".format(self.chunk_size))
        print("  > Max end datetime: {}".format(end_time_string))
        print("")
