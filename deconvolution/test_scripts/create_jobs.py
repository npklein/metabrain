#!/usr/bin/env python3

"""
File:         create_jobs.py
Created:      2020/04/22
Last Changed: 2020/04/26
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
import os

# Third party imports.

# Local application imports.

# Metadata
__program__ = "CellMap Profiles"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
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
        self.output_directory = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/"
        self.name_prefix = "CIA"
        self.time = "05:59:00"
        self.exclude = None
        self.cores = 1
        self.mem = 2
        self.count = 0
        self.start_index = 0
        self.stop_index = 21078
        self.batch_size = 50
        self.samples = 3703

    def start(self):
        start_indices = [i for i in range(self.start_index, self.stop_index, self.batch_size)]
        for job_id, start_index in enumerate(start_indices):
            batch_size = self.batch_size
            if job_id == (len(start_indices) - 1):
                batch_size = self.stop_index - start_index
            self.write_job_file(job_id, start_index, batch_size)

    def write_job_file(self, job_id, start_index, batch_size):
        skip_rows = ""
        if start_index > 0:
            skip_rows = " -sr {}".format(start_index)

        job_name = "{}{}".format(self.name_prefix, self.count + job_id)
        out_filepath = os.path.join(self.output_directory, 'output', job_name + ".out")
        bash_filepath = os.path.join(self.output_directory, job_name + ".sh")

        cores = ""
        if (self.cores - 1) > 0:
            cores = " -c {}".format(self.cores - 1)

        exclude = ""
        if self.exclude is not None:
            exclude = "#SBATCH --exclude={}"

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
                 "python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser.py{} -ne {} -ns {}{}\n".format(skip_rows, batch_size, self.samples, cores),
                 "\n",
                 "deactivate\n"]

        with open(bash_filepath, "w") as f:
            for line in lines:
                f.write(line)
        f.close()


if __name__ == '__main__':
    m = main()
    m.start()
