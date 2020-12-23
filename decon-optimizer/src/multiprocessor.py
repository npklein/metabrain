"""
File:         multiprocessor.py
Created:      2020/12/22
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
import multiprocessing as mp
import queue
import random
import time

# Third party imports.

# Local application imports.


class MultiProcessor:
    def __init__(self, samples, cores, max_end_time, log):
        self.samples = samples
        self.cores = cores
        self.max_end_time = max_end_time
        self.log = log

    def process(self, work_func, samples=None):
        if samples is None:
            samples = self.samples

        cores = self.cores
        if cores > len(samples):
            cores = len(samples)

        if cores == 1:
            result = self.singe_process(work_func=work_func,
                                        samples=samples)
        else:
            result = self.multiprocessing(work_func=work_func,
                                          samples=samples,
                                          cores=cores)

        missing = set(samples).symmetric_difference(set(result.keys()))
        if len(missing) > 0:
            self.log.info("\t\tMissing {} samples: '{}'.".format(len(missing), ",".join(missing)))

        return result

    def singe_process(self, work_func, samples):
        self.log.info("\t\tProcessing {} samples with 1 core.".format(len(samples)))

        result = {}
        start_time = int(time.time())
        for i, sample in enumerate(samples):
            if self.max_end_time <= time.time():
                break

            result[sample] = work_func(sample)

            rt_min, rt_sec = divmod(time.time() - start_time, 60)
            rt_hour, rt_min = divmod(rt_min, 60)
            self.log.info("\t\t\t[{:02d}:{:02d}:{:02d}] {}/{} samples completed [{:.2f}%]".format(int(rt_hour), int(rt_min), int(rt_sec), i+1, len(samples), (100 / len(samples)) * (i+1)))

        return result

    def multiprocessing(self, work_func, samples, cores):
        # Create queues.
        in_manager = mp.Manager()
        input_q = in_manager.Queue(maxsize=len(samples))
        out_manager = mp.Manager()
        output_q = out_manager.Queue(maxsize=len(samples))

        # Populate queues.
        for sample in samples:
            input_q.put(sample)

        # Create workers.
        processes = []
        for proc_id in range(cores):
            processes.append(mp.Process(target=self.worker,
                                        args=(input_q, output_q, work_func,
                                              self.max_end_time)
                                        ))

        # Start workers.
        started_procs = []
        for proc in processes:
            try:
                proc.start()
                started_procs.append(proc)
            except RuntimeError:
                pass
        del processes

        # Log the progress
        self.log.info("\t\tProcessing {} samples with {} cores.".format(len(samples), cores))
        n = None
        previous_n = n
        start_time = int(time.time())
        last_print_time = int(time.time())
        while n < len(samples):
            time.sleep(2)
            n = output_q.qsize()
            if self.max_end_time <= time.time():
                self.log.error("\t\t\tmax progress time reached")
                break

            now_time = int(time.time())
            if n != previous_n or now_time - last_print_time >= 30:
                last_print_time = now_time
                rt_min, rt_sec = divmod(time.time() - start_time, 60)
                rt_hour, rt_min = divmod(rt_min, 60)
                self.log.info("\t\t\t[{:02d}:{:02d}:{:02d}] {}/{} samples completed [{:.2f}%]".format(int(rt_hour), int(rt_min), int(rt_sec), n, len(samples), (100 / len(samples)) * n))
                previous_n = n

        # Stop the cores.
        for _ in range(cores):
            input_q.put("DONE")

        # Join the workers.
        for proc in started_procs:
            proc.join()

        # Join the results.
        result = {}
        while output_q.qsize() > 0:
            (sample, output) = output_q.get(True, timeout=1)
            result[sample] = output

        del in_manager, input_q, out_manager, output_q

        return result

    @staticmethod
    def worker(input_queue, output_queue, work_func, max_end_time):
        while True:
            if max_end_time <= time.time():
                break
            try:
                if not input_queue.empty():
                    sample = input_queue.get(True, timeout=1)
                    if sample == "DONE":
                        break

                    output_queue.put((sample, work_func(sample)))
            except queue.Empty:
                pass

            time.sleep(random.uniform(0.5, 1))
