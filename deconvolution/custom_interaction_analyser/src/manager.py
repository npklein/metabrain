"""
File:         manager.py
Created:      2020/04/01
Last Changed: 2020/04/10
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
import glob
import os

# Third party imports.

# Local application imports.
from general.local_settings import LocalSettings
from general.utilities import check_file_exists, prepare_output_dir
from .workers import process_worker


class Manager:
    """
    Class for the manager.
    """
    def __init__(self, settings_file, skip_rows, n_eqtls, n_samples, cores,
                 include, verbose):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param skip_rows: int, the number of rows to skip.
        :param n_eqtls: int, the number of eqtls in the input files.
        :param n_samples: int, the number of samples in the input files.
        :param cores: int, the number of cores to use.
        :param verbose: boolean, whether or not to print all update info.
        :param include: boolean, whether or not to include the unfinished
                        wait_list.
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
        self.analysis_max_runtime = settings.get_setting(
            "single_eqtl_max_runtime_in_min") * 60
        self.max_end_time = int(time.time()) + settings.get_setting(
            "max_runtime_in_hours") * 60 * 60
        self.panic_time = self.max_end_time - (settings.get_setting("panic_time_in_min") * 60)
        self.skip_rows = skip_rows
        self.verbose = verbose

        # Define output filenames.
        self.wait_list_outfile = os.path.join(self.outdir, settings.get_setting("wait_list_pickle_filename") + "{}.pkl".format(int(time.time())))
        self.perm_orders_outfile = os.path.join(self.outdir, settings.get_setting("permutations_order_pickle_filename") + ".pkl")
        self.pvalues_outfile = os.path.join(self.outdir, settings.get_setting("actual_pvalues_pickle_filename") + "{}.pkl".format(int(time.time())))
        self.perm_pvalues_outfile = os.path.join(self.outdir, settings.get_setting("permuted_pvalues_pickle_filename") + "{}.pkl".format(int(time.time())))

        # Load the previous wait list if need be.
        self.prev_wait_list = []
        if include:
            self.prev_wait_list = self.load_previous_wait_list(self.outdir, settings.get_setting("wait_list_pickle_filename"))

        # Count the number of eQTLs / samples if it isn't given by the user.
        if n_eqtls is None:
            n_eqtls = self.count_n_eqtls()
        self.n_eqtls = n_eqtls
        if n_samples is None:
            n_samples = self.count_n_samples()
        self.n_samples = n_samples

        # Determine optimal distribution of the work of the cores.
        self.n_cores, self.chunk_size = self.divide_work(n_eqtls, cores)

    def load_previous_wait_list(self, indir, filename):
        """
        Method that loads the wait list of a previous run.

        :param indir: string, the input directory containing the pickle files.
        :param filename: string, the prefix name of the input file.
        :return wait_list: list, the previous wait_lists.
        """
        # Order the filenames based on the integers appendices.
        infiles = glob.glob(os.path.join(indir, filename + "*.pkl"))
        start_indices = []
        for fpath in infiles:
            fpath_si = int(fpath.split(filename)[1].split('.')[0])
            start_indices.append(fpath_si)
        start_indices = sorted(start_indices)

        # Combine the found files.
        wait_list = []
        for i, start_index in enumerate(start_indices):
            fpath = os.path.join(indir, filename + str(start_index) + ".pkl")
            wait_list.extend(self.load_pickle(fpath))

        return wait_list

    def count_n_eqtls(self):
        """
        Method for counting the number of eQTLs (rows) in the genotype and
        expression file.

        :return geno_nrows: int, the number of eQTLs in both files.
        """
        geno_nrows = self.get_number_of_rows(self.geno_inpath)
        expr_nrows = self.get_number_of_rows(self.expr_inpath)
        if geno_nrows != expr_nrows:
            print("Number of rows in genotype / expression file do not match.")
            exit()

        return geno_nrows

    @staticmethod
    def get_number_of_rows(inpath):
        """
        Method for counting the number of rows in a .txt.gz file.

        :param inpath: string, the input file.
        :return i: int, the number of rows in the input file.
        """
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                pass
        f.close()
        return i

    def count_n_samples(self):
        """
        Method for counting the number of samples (columns) in the genotype and
        expression file.

        :return geno_ncols: int, the number of samples in both files.
        """
        geno_ncols = self.get_number_of_columns(self.geno_inpath)
        expr_ncols = self.get_number_of_columns(self.expr_inpath)
        if geno_ncols != expr_ncols:
            print("Number of columns in genotype / expression file do "
                  "not match.")
            exit()

        return geno_ncols

    @staticmethod
    def get_number_of_columns(inpath):
        """
        Method for counting the number of columns in a .txt.gz file.

        :param inpath: string, the input file.
        :return i: int, the number of columns in the input file.
        """
        with gzip.open(inpath, 'rb') as f:
            for i, l in enumerate(f):
                ncols = len(l.decode().split("\t")) - 1
                break
        f.close()
        return ncols

    @staticmethod
    def divide_work(n_eqtls, cores):
        """
        Method for distributing the work over the number of cores.
        .
        :param n_eqtls: int, the number of eqtls in the input files.
        :param cores: int, the number of cores to use.
        :return cores: int, the number of cores to use.
        :return chunk_size: int, the number of eQTLs per job entree.
        """
        # Cap the number of cores to the maximum.
        host_cores = mp.cpu_count()
        if cores > host_cores:
            cores = host_cores

        # Set the optimal values for the chunk_size and number of cores.
        if n_eqtls == 1:
            chunk_size = 1
            cores = 1
        elif (n_eqtls > 1) and math.floor(n_eqtls / cores) == 0:
            chunk_size = 1
            cores = n_eqtls
        else:
            chunk_size = max(1, min(math.ceil(n_eqtls / cores / 2), 10))

        return cores, chunk_size

    def start(self):
        """
        Method to start the manager.
        """
        self.print_arguments()
        print("[manager]\tstarting manager [{}]".format(
            datetime.now().strftime("%H:%M:%S")), flush=True)

        # Start the timer.
        start_time = int(time.time())

        # Get the permutation orders.
        permutation_orders = None
        if check_file_exists(self.perm_orders_outfile):
            print("[manager]\tloading permutation order.", flush=True)
            permutation_orders = self.load_pickle(self.perm_orders_outfile)

            # Validate the permutation orders for the given input.
            if len(permutation_orders) != (self.n_permutations + 1):
                print("[manager]\t\tinvalid.", flush=True)
                permutation_orders = None
            for order in permutation_orders:
                if len(order) != self.n_samples:
                    print("[manager]\t\tinvalid.", flush=True)
                    permutation_orders = None
                    break
        if permutation_orders is None:
            print("[manager]\tcreating permutation order.")
            permutation_orders = self.create_perm_orders()
            self.dump_pickle(permutation_orders, self.perm_orders_outfile)

        # Create the multiprocessing queues.
        start_manager = mp.Manager()
        start_q = start_manager.Queue()
        job_manager = mp.Manager()
        job_q = job_manager.Queue()
        result_manager = mp.Manager()
        result_q = result_manager.Queue()

        # Load the start queue.
        print("[manager]\tloading start queue.", flush=True)
        for _ in range(self.n_cores):
            start_q.put(permutation_orders)
        all_sample_orders = [i for i in range(len(permutation_orders))]

        # Load the wait list.
        print("[manager]\tcreating wait list.", flush=True)
        wait_list = self.load_wait_list(all_sample_orders)
        if self.prev_wait_list:
            print("[manager]\tadded {} eQTLs from previous "
                  "wait list.".format(len(self.prev_wait_list)), flush=True)
            wait_list.extend(self.prev_wait_list)

        # Start the workers.
        print("[manager]\tstarting processes.", flush=True)
        processes = []
        for worker_id in range(self.n_cores):
            processes.append(mp.Process(target=process_worker,
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
                                              self.verbose)))
        for proc in processes:
            proc.start()

        # Initialize variables.
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
        panic = False

        print("[manager]\tstarting scheduler.")
        while True:
            # Break the loop of the maximum runtime has been reached.
            # This prevents the workers to continue indefinitely.
            if time.time() > self.max_end_time:
                break

            # DOCTOR: check if a client isn't responding. If this is the case
            # then reschedule the unfinished work.
            now_time = int(time.time())
            tmp_doctor_dict = doctor_dict.copy()
            for worker_id, last_hr in tmp_doctor_dict.items():
                if (now_time - last_hr) > self.analysis_max_runtime:
                    # A dead worker has been found.
                    print("[doctor]\toh no, 'worker {}' "
                          "has died.".format(worker_id))
                    if worker_id not in dead_workers and \
                            worker_id in schedule.keys():
                        unfinished_work = schedule[worker_id]
                        if unfinished_work is not None:
                            # Reassign the unfinished work.
                            for key, value in unfinished_work.items():
                                wait_list.append((key, value, 1))

                            # Stop tracking this worker.
                            del schedule[worker_id]

                    if worker_id in doctor_dict.keys():
                        del doctor_dict[worker_id]
                    if worker_id not in dead_workers:
                        dead_workers.append(worker_id)

            # RECEIVER: Install new workers and reduce results.

            # Empty the results queue.
            while not result_q.empty():
                # Get an result from the queue.
                (worker_id, category, eqtl_index, sample_order_id,
                 data) = result_q.get()

                # Check the worker is not a zombie.
                if worker_id in dead_workers:
                    continue

                # Safe the moment the message was received. Don't save it
                # if it was already declared dead.
                if worker_id in doctor_dict.keys():
                    doctor_dict[worker_id] = int(time.time())

                # Check what kind of message it is.
                if category == "online":
                    print("[receiver]\ta new worker connected: "
                          "'worker {}'".format(worker_id),
                          flush=True)

                    # Add the worker to the schedule so it can receive work.
                    schedule[worker_id] = None
                    if not columns_added:
                        # Safe the columns.
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

                    # Prevent duplicates by checking if the worker didn't die
                    # already.
                    if worker_id in schedule.keys():
                        # Increment the counter.
                        counter += 1

                        # Safe the output to the right result list.
                        if sample_order_id == 0:
                            pvalue_data.append([eqtl_index] + data)
                        else:
                            perm_pvalues.extend(data[1:])

                        # Remove the job from the schedule.
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

            # Check if there is work to schedule.
            if len(wait_list) > 0:
                # Count the number of free_processes
                free_workers = []
                for worker_id, work in schedule.items():
                    if work is None:
                        free_workers.append(worker_id)

                # Check if there are free processes.
                if len(free_workers) > 0:
                    if len(wait_list) < len(free_workers):
                        free_workers = free_workers[0:len(wait_list)]

                    for worker_id in free_workers:
                        if self.verbose:
                            print("[scheduler]\tassigning work to "
                                  "'worker {}'".format(worker_id), flush=True)

                        # Assign work.
                        (start_index, sample_orders, chunk_size) = wait_list.pop(0)
                        schedule[worker_id] = {i: sample_orders.copy()
                                               for i in range(start_index,
                                                              start_index + chunk_size)}

            # Clear the job queue and push the new schedule. Push some
            # extra copies to make sure each process can get one.
            # clearing the queue prevents from clients executing outdated
            # commands.
            while not job_q.empty():
                job_q.get()
            extra_schedules = math.ceil(self.n_cores / 5)
            for _ in range(self.n_cores + extra_schedules):
                job_q.put(schedule)

            # END

            # Print statements to update the end-user.
            now_time = int(time.time())
            if now_time - last_print_time >= self.print_interval:
                last_print_time = now_time
                if self.verbose:
                    self.print_status(wait_list, schedule, eqtl_id_len, order_id_len)
                else:
                    self.print_progress(counter, total_analyses)

            # Check if we are done.
            if not wait_list:
                unfinished_work = schedule.values()
                unfinished_work = [x for x in unfinished_work if x]
                if not unfinished_work:
                    break

            # Check whether we are almost running out of time.
            if time.time() > self.panic_time:
                panic = True
                break

            time.sleep(self.sleep_time)

        # Update the end-user again. This allows for conformation all work
        # is indeed done.
        if self.verbose or panic:
            self.print_status(wait_list, schedule, eqtl_id_len, order_id_len)
        else:
            self.print_progress(counter, total_analyses)
        print("[manager]\tscheduler finished.", flush=True)

        # Send the kill signal.
        print("[manager]\tkilling processes.", flush=True)
        for proc in processes:
            proc.terminate()

        # Join the processes.
        for proc in processes:
            proc.join()

        # Pickle the output lists.
        print("[manager]\tsaving output lists.", flush=True)
        self.dump_pickle(pvalue_data, self.pvalues_outfile)
        self.dump_pickle(perm_pvalues, self.perm_pvalues_outfile)

        # Empty the schedule.
        print("[manager]\tclearing schedule, saving unfinished work.",
              flush=True)
        for unfinished_work in schedule.values():
            if unfinished_work is not None:
                for key, value in unfinished_work.items():
                    wait_list.append((key, value, 1))
        if wait_list:
            self.dump_pickle(wait_list, self.wait_list_outfile)

        # Print the process time.
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

        # Shutdown the manager.
        print("[manager]\tshutting down manager [{}]".format(
            datetime.now().strftime("%H:%M:%S")),
            flush=True)

    @staticmethod
    def load_pickle(fpath):
        """
        Method for loading data from a pickle file.

        :param fpath: string, the pickle input file.
        :return content: , the pickle content.
        """
        with open(fpath, "rb") as f:
            content = pickle.load(f)
        f.close()

        return content

    @staticmethod
    def dump_pickle(content, fpath):
        """
        Method for dumping data to a pickle file.

        :param content: , the pickle content.
        :param fpath: string, the pickle output file.
        """
        with open(fpath, "wb") as f:
            pickle.dump(content, f)
        f.close()

    def create_perm_orders(self):
        """
        Method for creating x random shuffles of the column indices.

        return sample_orders: list, a nested list with sample orders.
        """
        sample_orders = [self.get_column_indices()]
        for _ in range(self.n_permutations):
            sample_orders.append(self.create_random_shuffle())
        return sample_orders

    def get_column_indices(self):
        """
        Method for getting a ordered list of sample indices.

        :return : list, ordered list with sample indices.
        """
        return [x for x in range(self.n_samples)]

    def create_random_shuffle(self):
        """
        Method for creating a random shuffle of the column indices.

        :return column_indices: list, a shuffles list with sample indices.
        """
        column_indices = self.get_column_indices()
        random.shuffle(column_indices)
        return column_indices

    def load_wait_list(self, all_sample_orders):
        """
        Method for loading the wait list.

        :param all_sample_orders: list, a nested list with sample orders.
        :return wait_list: list, the wait list for the manager to work on.
        """
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
        """
        Method for printing the schedule in the terminal.

        :param wait_list: list, the wait list for the manager to work on.
        :param schedule: dict, the schedule of the manager.
        :param eqtl_id_len: int, the length of the highest eQTL ID.
        :param order_id_len: int, the length of the highest order ID.
        """
        print("[manager]", flush=True)
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
        print("[manager]", flush=True)

    @staticmethod
    def print_progress(counter, total_analyses):
        """
        Method for printing a progress counter in the terminal.

        :param counter: int, the number of analyses done.
        :param total_analyses: int, the total number of analyses to perform.
        """
        print_len = len(str(total_analyses))
        print("[manager]\treceived data from {:<{}d}/{:<{}d} analyses "
              "[{:.2f}%]".format(counter, print_len,
                                 total_analyses, print_len,
                                 (100 / total_analyses) * counter),
              flush=True)

    def print_arguments(self):
        """
        Method for printing the variables of the class.
        """
        end_time_string = datetime.fromtimestamp(self.max_end_time).strftime(
            "%d-%m-%Y, %H:%M:%S")
        panic_time_string = datetime.fromtimestamp(self.panic_time).strftime(
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
        print("  > Panic datetime: {}".format(panic_time_string))
        print("  > Skip rows: {}".format(self.skip_rows))
        print("  > EQTLs: {}".format(self.n_eqtls))
        print("  > Samples: {}".format(self.n_samples))
        print("  > Cores: {}".format(self.n_cores))
        print("  > Chunk size: {}".format(self.chunk_size))
        print("  > Verbose: {}".format(self.verbose))
        print("", flush=True)
