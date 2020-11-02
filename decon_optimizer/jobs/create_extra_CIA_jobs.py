#!/usr/bin/env python3

"""
File:         create_extra_CIA_jobs.py
Created:      2020/10/25
Last Changed: 2020/10/27
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
import argparse
import os

# Third party imports.

# Local application imports.

# Metadata
__program__ = "Create Extra CIA jobs"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.job = getattr(arguments, 'job').upper()
        self.input_folder = getattr(arguments, 'input').lower()
        self.output_folder = getattr(arguments, 'output').lower()
        self.settings = getattr(arguments, 'settings').lower()
        self.exclude = getattr(arguments, 'exclude')
        self.batch_size = getattr(arguments, 'batch')
        self.n_samples = getattr(arguments, 'n_samples')
        self.cores = 1
        self.mem = getattr(arguments, 'mem')

        # Hard coded.
        self.missing = ['18-19', '58-59', '117-119', '145-149', '155-159', '235-239', '245-249', '254-259', '335-339', '366-369', '376-379', '1727-1729', '1739', '1749', '1758-1759', '1769', '1779', '1789', '1799', '1809', '1819', '1829', '1839', '1849', '1859', '1869', '1878-1879', '1889', '1899', '1909', '1919', '1929', '1939-1949', '1959', '1969', '1979', '1989', '1999', '2638-2639', '2649', '4538-4539', '4549', '4559', '4569', '4577-4579', '4588-4589', '4598-4609', '4619', '4629', '4639', '4649', '4659-4669', '4679', '4689', '4699', '4709', '4719', '4729', '4739', '4749', '4759', '4769', '4779', '4789', '4799', '4809-4819', '4827-4829', '4837-4839', '4847-4849', '4857-4859', '4867-4869', '4877-4879', '4887-4889', '4907-4909', '4917-4919', '4927-4929', '4967-4969', '4977-4979', '5017-5019', '5047-5049', '5117-5119', '5157-5159', '5170-5179', '5207-5209', '5227-5229', '5247-5249', '5257-5259', '5267-5269', '5277-5279', '5297-5299', '5316-5319', '5327-5329', '5338-5339', '5376-5379', '5396-5399', '5416-5419', '5426-5429', '5436-5439', '5446-5449', '5456-5459', '5466-5469', '5476-5479', '5486-5489', '5496-5499', '5506-5509', '5516-5519', '5526-5529', '5536-5539', '5546-5549', '5556-5559', '5566-5569', '5576-5579', '5586-5589', '5596-5599', '5616-5619', '5626-5629', '5636-5639', '5646-5649', '5656-5659', '5666-5669', '5676-5679', '5686-5689', '5697-5709', '6397-6399', '6458-6459', '6959', '8219-8229', '8238-8239', '8248-8249', '8258-8259', '8268-8269', '8278-8279', '8288-8289', '8298-8299', '8308-8309', '8318-8319', '8328-8329', '8338-8339', '8348-8349', '8358-8359', '8368-8369', '8378-8379', '8388-8389', '8398-8399', '8418-8419', '8428-8429', '8438-8439', '8458-8459', '8469', '8508-8509', '8528-8529', '8548-8549', '8807-8809', '8889', '9719', '10060-10069', '10209', '11220-11229', '11250-11259']

        # Set the variables.
        self.outdir = os.path.join(Path(__file__).parent.absolute(), self.job + "_E")
        self.log_file_outdir = os.path.join(self.outdir, 'output')
        time = getattr(arguments, 'time').lower()
        self.time = None
        if time == "short":
            self.time = "05:59:00"
        elif time == "medium":
            self.time = "23:59:00"
        elif time == "long":
            self.time = "6-23:49:00"

        for outdir in [self.outdir, self.log_file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-j",
                            "--job",
                            type=str,
                            required=True,
                            help="The name of the job.")
        parser.add_argument("-i",
                            "--input",
                            type=str,
                            default=None,
                            help="The name of the input directory. "
                                 "Default: 'output'.")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            default=None,
                            help="The name of the output directory. "
                                 " Default: same as -i / --input.")
        parser.add_argument("-s",
                            "--settings",
                            type=str,
                            required=True,
                            help="The settings input file (without '.json'), "
                                 "default: 'default_settings'.")
        parser.add_argument("-b",
                            "--batch",
                            type=int,
                            default=50,
                            help="The number of eQTLs per job. Default: 50.")
        parser.add_argument("-ns",
                            "--n_samples",
                            type=int,
                            required=True,
                            help="The number of samples.")
        parser.add_argument("-t",
                            "--time",
                            type=str,
                            default="short",
                            choices=["short", "medium", "long"],
                            help="The time required.")
        parser.add_argument("-c",
                            "--cores",
                            type=int,
                            default=1,
                            help="The number of cores required.")
        parser.add_argument("-m",
                            "--mem",
                            type=int,
                            default=2,
                            help="The memory required.")
        parser.add_argument("-e",
                            "--exclude",
                            type=str,
                            default=None,
                            help="The name of node to exclude,"
                                 "default: None.")


        return parser.parse_args()

    def start(self):
        jobs_info = self.create_jobs_info()

        job_id_count = 0
        for (start_index, batch_size) in jobs_info:
            self.write_job_file(job_id_count, start_index, batch_size)
            job_id_count += 1

    def create_jobs_info(self):
        jobs_info = []
        for missing in self.missing:
            split_missing = missing.split("-")

            start_index = int(split_missing[0])
            if len(split_missing) == 1:
                stop_index = start_index + 1
            elif len(split_missing) == 2:
                stop_index = int(split_missing[1]) + 1
            else:
                print("unexpected input")
                continue

            batch_size = stop_index - start_index

            if batch_size > self.batch_size:
                subjobs = [(i, self.batch_size) for i in range(start_index,
                                                               stop_index,
                                                               self.batch_size)]
                if (stop_index - subjobs[-1][0]) != self.batch_size:
                    subjobs[-1] = (subjobs[-1][0], stop_index - subjobs[-1][0])
            else:
                subjobs = [(start_index, batch_size)]

            jobs_info.extend(subjobs)

        print(jobs_info)

        return jobs_info

    def write_job_file(self, job_id, start_index, batch_size):
        skip_rows = ""
        if start_index > 0:
            skip_rows = " -sr {}".format(start_index)

        job_name = "{}_E{}".format(self.job.upper(), job_id)
        out_filepath = os.path.join(self.log_file_outdir, job_name + ".out")
        bash_filepath = os.path.join(self.outdir, job_name + ".sh")

        exclude = ""
        if self.exclude is not None:
            exclude = "#SBATCH --exclude={}".format(self.exclude)

        lines = ["#!/bin/bash\n",
                 "#SBATCH --job-name={}\n".format(job_name),
                 "#SBATCH --output={}\n".format(out_filepath),
                 "#SBATCH --error={}\n".format(out_filepath),
                 "#SBATCH --time={}\n".format(self.time),
                 "#SBATCH --cpus-per-task={}\n".format(self.cores),
                 "#SBATCH --mem={}gb\n".format(self.mem),
                 "#SBATCH --nodes=1\n",
                 "#SBATCH --open-mode=append\n",
                 "#SBATCH --export=NONE\n",
                 "#SBATCH --get-user-env=L\n",
                 "{}\n".format(exclude),
                 "\n",
                 "module load Python/3.6.3-foss-2015b\n",
                 "source $HOME/venv/bin/activate\n",
                 "\n",
                 "python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-decon-optimizer/custom_interaction_analyser.py -i {} -o {} -s {}{} -ne {} -ns {}\n".format(self.input_folder, self.output_folder, self.settings, skip_rows, batch_size, self.n_samples),
                 "\n",
                 "deactivate\n"]

        with open(bash_filepath, "w") as f:
            for line in lines:
                f.write(line)
        f.close()
        print("Created {}".format(bash_filepath))


if __name__ == '__main__':
    m = main()
    m.start()
