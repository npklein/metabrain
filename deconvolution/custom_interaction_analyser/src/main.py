"""
File:         main.py
Created:      2020/03/23
Last Changed: 2020/03/26
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
import queue
import gzip
import time
import math
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
from general.utilities import prepare_output_dir
from general.df_utilities import save_dataframe
from .workers import process_worker


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, n_eqtls, cores):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param n_eqtls: int, the number of eqtls in the input file.
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
        self.max_end_time = int(time.time()) + settings.get_setting("max_runtime_in_min") * 60
        self.n_permutations = settings.get_setting("n_permutations")

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Count the number of eQTLs.
        if n_eqtls is None:
            n_eqtls = self.count_n_eqtls()

        # Cap the number of cores to the maximum.
        host_cores = mp.cpu_count()
        if cores > host_cores:
            cores = host_cores

        # Set the correct values for chunk size, and cores.
        if n_eqtls == 1:
            chunk_size = 1
            cores = 1
        elif (n_eqtls > 1) and math.floor(n_eqtls / cores) == 0:
            chunk_size = 1
            cores = n_eqtls
        else:
            chunk_size = math.ceil(n_eqtls / cores / 2)
        self.n_cores = cores
        self.n_eqtls = n_eqtls
        self.chunk_size = chunk_size

    def count_n_eqtls(self):
        print("Finding the number of eQTLs")
        geno_nrows = self.get_number_of_lines(self.geno_inpath)
        expr_nrows = self.get_number_of_lines(self.expr_inpath)
        if geno_nrows != expr_nrows:
            print("Number of rows in genotype / expression file do not match.")
            exit()

        return geno_nrows

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()

        # Create queues.
        input_manager = mp.Manager()
        input_queue = input_manager.Queue()
        output_manager = mp.Manager()
        output_queue = output_manager.Queue()


        if self.n_eqtls == 1 and self.chunk_size == 1:
            input_queue.put(1)
        else:
            for i in range(1, self.n_eqtls + 1, self.chunk_size):
                input_queue.put(i)

        # Fire off the workers.
        processes = []
        for worker_id in range(self.n_cores):
            processes.append(mp.Process(target=process_worker,
                                        args=(worker_id,
                                              self.cov_inpath,
                                              self.geno_inpath,
                                              self.expr_inpath,
                                              self.tech_covs,
                                              input_queue,
                                              output_queue,
                                              self.max_end_time,
                                              self.chunk_size,
                                              self.n_permutations)))
        for proc in processes:
            proc.start()

        # Empty the output queue.
        columns = None
        eqtl_zscores = []
        eqtl_adj_zscores = []
        last_print_time = time.time()
        while (len(eqtl_zscores) != self.n_eqtls) or (len(eqtl_adj_zscores) != self.n_eqtls):
            now_time = time.time()
            if self.max_end_time <= now_time:
                print("max progress time reached")
                exit()
            if now_time - last_print_time >= 60:
                last_print_time = now_time
                print("\tReceived data from {}/{} eQTLs "
                      "[{:.2f}%]".format(len(eqtl_zscores),
                                         self.n_eqtls,
                                         (100 / self.n_eqtls) * len(eqtl_zscores)))

            try:
                if output_queue.empty():
                    time.sleep(0.1)
                    continue

                # Get the new output
                (zscore_type, new_output) = output_queue.get(True, timeout=1)
                new_index = new_output[0]

                # Check the length.
                if columns is not None and len(new_output) != len(columns):
                    print("Length of new ouput on index '{}' does not "
                          "correspond with the remainder of the "
                          "data.".format(new_index))
                    exit()

                # Safe the columns.
                if new_index == -1:
                    if columns is None:
                        columns = new_output
                    else:
                        if new_output != columns:
                            print("Columns are not identical over "
                                  "subprocesses.")
                else:
                    # Safe the zscores in any position.
                    if zscore_type == "default":
                        eqtl_zscores.append(new_output)
                    elif zscore_type == "adjusted":
                        eqtl_adj_zscores.append(new_output)

            except queue.Empty:
                time.sleep(1)
                continue

        # Join the processes.
        for proc in processes:
            proc.join()

        # Sort the dataframes.
        inter_df = self.create_df(eqtl_zscores, columns)
        adj_inter_df = self.create_df(eqtl_adj_zscores, columns)

        # Pickle the files.
        inter_df.to_pickle(os.path.join(self.outdir, "inter_df.pkl"))
        adj_inter_df.to_pickle(os.path.join(self.outdir, "adj_inter_df.pkl"))

        # inter_df = pd.read_pickle(os.path.join(self.outdir, "inter_df.pkl"))
        # adj_inter_df = pd.read_pickle(os.path.join(self.outdir, "adj_inter_df.pkl"))

        print(inter_df)
        print(adj_inter_df)

        # Write the output file.
        save_dataframe(df=adj_inter_df,
                       outpath=os.path.join(self.outdir,
                                            "interaction_table.txt.gz"),
                       header=True, index=True)

        # Compare the two values
        self.compare_z_scores(inter_df, adj_inter_df, self.outdir)

        # Stop timer.
        process_time = time.time() - start_time

        # Print the time.
        minutes, seconds = divmod(process_time, 60)
        print("Finished in {:.0f} minutes and "
              "{:.0f} seconds".format(minutes, seconds))

    @staticmethod
    def get_number_of_lines(inpath):
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                pass
        return i

    @staticmethod
    def create_df(data, columns):
        df = pd.DataFrame(data, columns=columns)
        df.set_index(df.columns[0], inplace=True)
        df.sort_index(inplace=True)
        df.set_index("-", inplace=True)
        df = df.T

        return df

    @staticmethod
    def compare_z_scores(zscores_df, adj_zscores_df, outdir):
        df = pd.DataFrame({"default": zscores_df.melt()["value"],
                           "adjusted": adj_zscores_df.melt()["value"]})
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.dropna(inplace=True)
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
        g.set_ylabel('Adjusted Z-score',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('Default Z-score',
                     fontsize=8,
                     fontweight='bold')
        g.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        g.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        fig.savefig(os.path.join(outdir, "zscore_comparison.png"))

    def print_arguments(self):
        end_time_string = datetime.fromtimestamp(self.max_end_time).strftime("%d-%m-%Y, %H:%M:%S")
        print("Arguments:")
        print("  > Genotype datafile: {}".format(self.geno_inpath))
        print("  > Expression datafile: {}".format(self.expr_inpath))
        print("  > Covariates datafile: {}".format(self.cov_inpath))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > N eQTLs: {}".format(self.n_eqtls))
        print("  > N cores: {}".format(self.n_cores))
        print("  > Chunk size: {}".format(self.chunk_size))
        print("  > Max end datetime: {}".format(end_time_string))
        print("")
