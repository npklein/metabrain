"""
File:         manager.py
Created:      2020/04/01
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
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
import pickle
import random
import gzip
import time
import math
import os

# Third party imports.

# Local application imports.
from general.local_settings import LocalSettings
from general.utilities import prepare_output_dir
from .workers2 import process_worker


class Manager:
    """
    Class for the manager.
    """

    def __init__(self, settings_file, skip_rows, n_eqtls, n_samples, cores,
                 verbose):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param skip_rows: int, the number of rows to skip.
        :param n_eqtls: int, the number of eqtls in the input files.
        :param n_samples: int, the number of samples in the input files.
        :param cores: int, the number of cores to use.
        :param verbose: boolean, whether or not to print all update info.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Safe settings.
        self.geno_inpath = settings.get_setting("genotye_datafile")
        self.expr_inpath = settings.get_setting("expression_datafile")
        self.cov_inpath = settings.get_setting("covariates_datafile")
        self.tech_covs = settings.get_setting("technical_covariates")
        self.n_permutations = settings.get_setting("n_permutations")
        self.sleep_time = settings.get_setting("sleep_time_in_sec")
        self.print_interval = settings.get_setting("print_interval_in_sec")
        self.analysis_max_runtime = settings.get_setting("single_eqtl_max_runtime_in_min") * 60
        self.max_end_time = int(time.time()) + settings.get_setting(
            "max_runtime_in_hours") * 60 * 60
        self.verbose = verbose

        # Count the number of eQTLs / samples.
        if n_eqtls is None:
            n_eqtls = self.count_n_eqtls()
        if n_samples is None:
            n_samples = self.count_n_samples()

        # Determine optimal distribution.
        self.skip_rows, self.n_eqtls, self.n_samples, self.n_cores, \
        self.chunk_size = self.determine_optimal_distribution(skip_rows,
                                                              n_eqtls,
                                                              n_samples,
                                                              cores)

        # Define output filenames.
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

    @staticmethod
    def determine_optimal_distribution(skip_rows, n_eqtls, n_samples, cores):
        # Cap the number of cores to the maximum.
        host_cores = mp.cpu_count()
        if cores > host_cores:
            cores = host_cores

        # Set the correct values for skip_rows, chunk size, and cores.
        if n_eqtls == 1:
            chunk_size = 1
            cores = 1
        elif (n_eqtls > 1) and math.floor(n_eqtls / cores) == 0:
            chunk_size = 1
            cores = n_eqtls
        else:
            chunk_size = min(math.ceil(n_eqtls / cores / 2), 10)
        return skip_rows, n_eqtls, n_samples, cores, chunk_size

    def start(self):
        """
        Function to start the manager.
        """
        self.print_arguments()
        print("[manager]\tstarting manager [{}]".format(
            datetime.now().strftime("%H:%M:%S")), flush=True)

        # start the timer.
        start_time = int(time.time())

        # create queues.
        start_manager = mp.Manager()
        start_q = start_manager.Queue()
        job_manager = mp.Manager()
        job_q = job_manager.Queue()
        result_manager = mp.Manager()
        result_q = result_manager.Queue()

        # load the start queue.
        print("[manager]\tloading start queue.", flush=True)
        all_sample_orders = self.fill_start_queue(start_q)

        # load the wait list.
        print("[manager]\tcreating waitlist.", flush=True)
        wait_list = self.load_wait_list(all_sample_orders)

        # start the workers.
        print("[manager]\tstarting processes.", flush=True)
        processes = {}
        for worker_id in range(self.n_cores):
            processes[worker_id] = mp.Process(target=process_worker,
                                              args=(worker_id,
                                                    self.cov_inpath,
                                                    self.geno_inpath,
                                                    self.expr_inpath,
                                                    self.tech_covs,
                                                    start_q,
                                                    job_q,
                                                    result_q,
                                                    self.sleep_time,
                                                    self.max_end_time,
                                                    self.verbose))
        for proc in processes.values():
            proc.start()

        # initialize variables.
        dead_workers = []
        schedule = {}
        doctor_dict = {}
        columns_added = False
        pvalue_data = []
        perm_pvalues = []
        eqtl_id_len = len(str(self.n_eqtls))
        order_id_len = len(str(len(all_sample_orders)))
        last_print_time = int(time.time()) - self.print_interval
        counter = 0
        total_analyses = self.n_eqtls * (self.n_permutations + 1)

        print("[manager]\tstarting scheduler.")
        while True:
            if time.time() > self.max_end_time:
                break

            # DOCTOR: checks heartbeats

            # check if a client isn't responding.

            now_time = int(time.time())
            tmp_doctor_dict = doctor_dict.copy()

            for worker_id, last_hr in tmp_doctor_dict.items():
                if (now_time - last_hr) > self.analysis_max_runtime:
                    # a dead client has been found.
                    print("[doctor]\toh no, 'worker {}' "
                          "has died.".format(worker_id))
                    if worker_id not in dead_workers and \
                            worker_id in schedule.keys():
                        unfinished_work = schedule[worker_id]
                        if unfinished_work is not None:
                            # reassign the unfinished work.
                            for key, value in unfinished_work.items():
                                wait_list.append((key, value, 1))

                            # Stop tracking this worker.
                            del schedule[worker_id]

                    if worker_id in doctor_dict.keys():
                        del doctor_dict[worker_id]

                    dead_workers.append(worker_id)

            # RECEIVER: install new clients and reduce results

            # empty the results queue.
            while not result_q.empty():
                # get an result from the queue.
                (worker_id, category, eqtl_index, sample_order_id, data) = result_q.get()

                # safe the moment the message was received.
                if worker_id not in dead_workers:
                    doctor_dict[worker_id] = int(time.time())

                # check what kind of message it is.
                if category == "online":
                    print("[receiver]\ta new worker connected: "
                          "'worker {}'".format(worker_id),
                          flush=True)

                    schedule[worker_id] = None
                    if not columns_added:
                        pvalue_data.append([eqtl_index] + data)
                        columns_added = True
                elif category == "result":
                    if self.verbose:
                        print("[receiver]\treceived an output "
                              "(eQTL_{:<{}d}, order_{:<{}d}) from: "
                              "'worker {}'".format(eqtl_index,
                                                   eqtl_id_len,
                                                   sample_order_id,
                                                   order_id_len,
                                                   worker_id),
                              flush=True)

                    # prevent duplicates by checking if the worker didn't die
                    # already.
                    if worker_id in schedule.keys():
                        # Increment the counter.
                        counter += 1

                        # Safe the output to the right result.
                        if sample_order_id == 0:
                            pvalue_data.append([eqtl_index] + data)
                        else:
                            perm_pvalues.extend(data[1:])

                        # remove the job from the schedule.
                        if eqtl_index in schedule[worker_id].keys():
                            sample_orders_to_perform = schedule[worker_id][eqtl_index]
                            sample_orders_to_perform.remove(sample_order_id)
                            if sample_orders_to_perform:
                                schedule[worker_id][eqtl_index] = sample_orders_to_perform
                            else:
                                del schedule[worker_id][eqtl_index]

                            # remove the job so new work can be added.
                            if not schedule[worker_id].keys():
                                schedule[worker_id] = None
                else:
                    pass

            # SCHEDULER: map jobs

            # check if there is work to schedule.
            if len(wait_list) > 0:

                # count the number of free_processes
                free_workers = []
                for worker_id, work in schedule.items():
                    if work is None:
                        free_workers.append(worker_id)

                # check if there are free processes.
                if len(free_workers) > 0:
                    if len(wait_list) < len(free_workers):
                        free_workers = free_workers[0:len(wait_list)]

                    for worker_id in free_workers:
                        (start_index, sample_orders, chunk_size) = wait_list.pop(0)
                        if self.verbose:
                            print("[scheduler]\tassigning work to "
                                  "'worker {}'".format(worker_id), flush=True)
                        schedule[worker_id] = {i: sample_orders.copy()
                                               for i in range(start_index,
                                                              start_index + chunk_size)}

            # clear the job queue and push the new schedule. Push some
            # extra copies to make sure each process can get one.
            # clearing the queue prevents from clients executing outdated
            # commands.
            while not job_q.empty():
                job_q.get()
            extra_schedules = math.ceil(self.n_cores / 5)
            for _ in range(self.n_cores + extra_schedules):
                job_q.put(schedule)

            # END

            # print statements to update to end-user.
            now_time = int(time.time())
            if now_time - last_print_time >= self.print_interval:
                last_print_time = now_time
                if self.verbose:
                    self.print_status(wait_list, schedule, eqtl_id_len, order_id_len)
                else:
                    self.print_progress(counter, total_analyses)

            # check if we are done.
            if not wait_list:
                unfinished_work = schedule.values()
                unfinished_work = [x for x in unfinished_work if x]
                if not unfinished_work:
                    break

            time.sleep(self.sleep_time)

        # print statements to update to end-user.
        if self.verbose:
            self.print_status(wait_list, schedule, eqtl_id_len, order_id_len)
        else:
            self.print_progress(counter, total_analyses)
        print("[manager]\tscheduler finished.", flush=True)

        # Send the kill signal.
        print("[manager]\tkilling processes.", flush=True)
        for proc in processes.values():
            proc.terminate()

        # Join the processes.
        for proc in processes.values():
            proc.join()

        # Pickle the files.
        print("[manager]\tsaving output files.", flush=True)
        with open(self.pvalues_outfile, "wb") as f:
            pickle.dump(pvalue_data, f)
        f.close()
        with open(self.perm_pvalues_outfile, "wb") as f:
            pickle.dump(perm_pvalues, f)
        f.close()

        # Print the time.
        run_time = int(time.time()) - start_time
        run_time_min, run_time_sec = divmod(run_time, 60)
        run_time_hour, run_time_min = divmod(run_time_min, 60)
        print("[manager]\tfinished in  {} hour(s), {} minute(s) and "
              "{} second(s).".format(int(run_time_hour),
                                     int(run_time_min),
                                     int(run_time_sec)),
              flush=True)
        print("[receiver]\treceived {:.2f} analyses "
              "per second.".format(counter / run_time),
              flush=True)

        # shutdown the manager.
        print("[manager]\tshutting down manager [{}]".format(
            datetime.now().strftime("%H:%M:%S")),
            flush=True)

    def fill_start_queue(self, start_queue):
        # Create x random shuffles of the column indices.
        sample_orders = [self.get_column_indices()]
        for _ in range(self.n_permutations):
            sample_orders.append(self.create_random_shuffle())

        for _ in range(self.n_cores):
            start_queue.put(sample_orders)

        return [i for i in range(len(sample_orders))]

    def get_column_indices(self):
        return [x for x in range(self.n_samples)]

    def create_random_shuffle(self):
        column_indices = self.get_column_indices()
        random.shuffle(column_indices)
        return column_indices

    def load_wait_list(self, all_sample_orders):
        wait_list = []
        if self.n_eqtls == 1 and self.chunk_size == 1:
            wait_list.append((self.skip_rows + 1,
                              all_sample_orders.copy(),
                              self.chunk_size))
        else:
            for start_index in range(self.skip_rows + 1,
                                     self.skip_rows + self.n_eqtls + 1,
                                     self.chunk_size):
                chunk_size = self.chunk_size
                if (start_index + self.chunk_size) > (self.skip_rows + self.n_eqtls):
                    chunk_size = self.skip_rows + self.n_eqtls - start_index + 1

                wait_list.append((start_index, all_sample_orders.copy(), chunk_size))

        return wait_list

    @staticmethod
    def print_status(wait_list, schedule, eqtl_id_len, order_id_len):
        print("[manager]\t", flush=True)
        if len(wait_list) > 0:
            work_info = []
            wait_list = sorted(wait_list, key=lambda x: x[0])
            for job in wait_list:
                work_info.append("eQTL_{:<{}d}-eQTL_{:<{}d} (x{:<{}d})".format(job[0], eqtl_id_len,  job[0] + job[2], eqtl_id_len, len(job[1]), order_id_len))
            print("[manager]\twait list: [{}]".format(', '.join(work_info)),
                  flush=True)
        else:
            print("[manager]\twait list: empty")
        print("[manager]\tcurrent schedule:")
        if len(schedule.keys()) > 0:
            connected_workers = list(schedule.keys())
            connected_workers.sort()
            for worked_id in connected_workers:
                work_info = []
                work = schedule[worked_id]
                if work:
                    for key in sorted(work.keys()):
                        work_info.append(
                            "eQTL_{:<{}d} (x{:<{}d})".format(key, eqtl_id_len, len(schedule[worked_id][key]), order_id_len))
                print("[manager]\t\tworker {}: [{}]".format(worked_id, ', '.join(work_info)),
                      flush=True)
        else:
            print("[manager]\t\tempty", flush=True)

    @staticmethod
    def print_progress(counter, total_analyses):
        print_len = len(str(total_analyses))
        print("[manager]\treceived data from {:<{}d}/{:<{}d} analyses "
              "[{:.2f}%]".format(counter, print_len,
                                 total_analyses, print_len,
                                 (100 / total_analyses) * counter),
              flush=True)

    def print_arguments(self):
        end_time_string = datetime.fromtimestamp(self.max_end_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        print("Arguments:")
        print("  > Genotype datafile: {}".format(self.geno_inpath))
        print("  > Expression datafile: {}".format(self.expr_inpath))
        print("  > Covariates datafile: {}".format(self.cov_inpath))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > Actual P-values output file: {}".format(self.pvalues_outfile))
        print("  > Permutated P-values output file: {}".format(self.perm_pvalues_outfile))
        print("  > Permutations: {}".format(self.n_permutations))
        print("  > Sleep time: {} sec".format(self.sleep_time))
        print("  > Print interval: {} sec".format(self.print_interval))
        print("  > Single analysis max runtime: {} sec".format(self.analysis_max_runtime))
        print("  > Max end datetime: {}".format(end_time_string))
        print("  > Skip rows: {}".format(self.skip_rows))
        print("  > EQTLs: {}".format(self.n_eqtls))
        print("  > Samples: {}".format(self.n_samples))
        print("  > Cores: {}".format(self.n_cores))
        print("  > Chunk size: {}".format(self.chunk_size))
        print("  > Verbose: {}".format(self.verbose))
        print("", flush=True)
