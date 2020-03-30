"""
File:         main.py
Created:      2020/03/23
Last Changed: 2020/03/30
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
import pickle
import random
import queue
import gzip
import time
import math
import os

# Third party imports.

# Local application imports.
from general.local_settings import LocalSettings
from general.utilities import prepare_output_dir
from .workers import process_worker


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, skip_rows, n_eqtls, n_samples, cores):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param skip_rows: int, the number of rows to skip.
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

        # Set the correct values for skip_rows, chunk size, and cores.
        if n_eqtls == 1:
            chunk_size = 1
            cores = 1
        elif (n_eqtls > 1) and math.floor(n_eqtls / cores) == 0:
            chunk_size = n_eqtls
            cores = n_eqtls
        else:
            chunk_size = math.ceil(n_eqtls / cores / 2)
        self.skip_rows = skip_rows
        self.n_cores = cores
        self.n_eqtls = n_eqtls
        self.n_samples = n_samples
        self.chunk_size = chunk_size

        # Create output filepaths.
        self.pvalues_outfile = os.path.join(self.outdir, settings.get_setting("actual_pvalues_pickle_filename") + "{}.pkl".format(self.skip_rows))
        self.perm_pvalues_outfile = os.path.join(self.outdir, settings.get_setting("permuted_pvalues_pickle_filename") + "{}.pkl".format(self.skip_rows))

    def count_n_eqtls(self):
        print("Finding the number of eQTLs")
        geno_nrows = self.get_number_of_rows(self.geno_inpath)
        expr_nrows = self.get_number_of_rows(self.expr_inpath)
        if geno_nrows != expr_nrows:
            print("Number of rows in genotype / expression file do not match.")
            exit()

        return geno_nrows

    @staticmethod
    def get_number_of_rows(inpath):
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                pass
        f.close()
        return i

    def count_n_samples(self):
        print("Finding the number of samples")
        geno_ncols = self.get_number_of_columns(self.geno_inpath)
        expr_ncols = self.get_number_of_columns(self.expr_inpath)
        if geno_ncols != expr_ncols:
            print("Number of columns in genotype / expression file do "
                  "not match.")
            exit()

        return geno_ncols

    @staticmethod
    def get_number_of_columns(inpath):
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                ncols = len(l.decode().split("\t")) - 1
                break
        f.close()
        return ncols

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

        # Fill input queue.
        print("Loading input queue.")
        njobs = self.fill_input_queue(input_queue)
        print("Number of jobs to complete: {}".format(njobs))

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
        columns_added = False
        pvalue_data = []
        perm_pvalues = []
        analyses_done = 0
        total_analyses = self.n_eqtls * (self.n_permutations + 1)
        print_len = len(str(total_analyses))
        last_print_time = time.time() - 50
        print("Waiting for output from workers.")
        while analyses_done != total_analyses:
            now_time = time.time()
            if self.max_end_time <= now_time:
                print("ERROR: max progress time reached")
                exit()
            if now_time - last_print_time >= 60:
                last_print_time = now_time
                print("\tReceived data from {:{}d}/{:{}d} analyses "
                      "[{:.2f}%]".format(analyses_done, print_len,
                                         total_analyses, print_len,
                                         (100 / total_analyses) * analyses_done))

            try:
                if output_queue.empty():
                    time.sleep(0.1)
                    continue

                # Get the new output
                (output_type, output) = output_queue.get(True, timeout=1)
                new_index = output[0]

                # Safe the columns.
                if new_index == -1 and not columns_added:
                    pvalue_data.append(output)
                    columns_added = True
                else:
                    # Safe the output to the right result.
                    analyses_done += 1
                    if output_type == -1:
                        pvalue_data.append(output)
                    else:
                        perm_pvalues.extend(output[2:])

            except queue.Empty:
                time.sleep(1)
                continue
        print("\tReceived data from {:{}d}/{:{}d} analyses "
              "[{:.2f}%]".format(analyses_done, print_len,
                                 total_analyses, print_len,
                                 (100 / total_analyses) * analyses_done))

        # Join the processes.
        for proc in processes:
            proc.join()

        # Pickle the files.
        with open(self.pvalues_outfile, "wb") as f:
            pickle.dump(pvalue_data, f)
        f.close()
        with open(self.perm_pvalues_outfile, "wb") as f:
            pickle.dump(perm_pvalues, f)
        f.close()

        # Stop timer.
        process_time = time.time() - start_time

        # Print the time.
        minutes, seconds = divmod(process_time, 60)
        print("Finished in {:.0f} minutes and "
              "{:.0f} seconds".format(minutes, seconds))

    def fill_input_queue(self, input_queue):
        # Create a list of start indices.
        start_indices = [i for i in range(self.skip_rows + 1,
                                          self.n_eqtls + 1,
                                          self.chunk_size)]
        if self.n_eqtls == 1 and self.chunk_size == 1:
            start_indices = [self.skip_rows + 1]

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

        return njobs

    def get_column_indices(self):
        return [x for x in range(self.n_samples)]

    def create_random_shuffle(self):
        column_indices = self.get_column_indices()
        random.shuffle(column_indices)
        return column_indices

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
        print("  > Actual P-values output file: {}".format(self.pvalues_outfile))
        print("  > Permutated P-values output file: {}".format(self.perm_pvalues_outfile))
        print("  > Skip rows: {}".format(self.skip_rows))
        print("  > N eQTLs: {}".format(self.n_eqtls))
        print("  > N samples: {}".format(self.n_samples))
        print("  > N cores: {}".format(self.n_cores))
        print("  > Chunk size: {}".format(self.chunk_size))
        print("  > Max end datetime: {}".format(end_time_string))
        print("")
