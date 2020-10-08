"""
File:         combine_eqtlprobes.py
Created:      2020/10/08
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
import os

# Third party imports.
import pandas as pd

# Local application imports.
from matrix_preparation.src.utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class CombineEQTLProbes:
    def __init__(self, settings, force, outdir):
        self.indir = settings["input_directory"]
        self.iter_dirname = settings["iteration_dirname"]
        self.n_iterations = settings["iterations"]
        self.in_filename = settings["in_filename"]
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'combine_eqtlprobes')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "eQTLprobes_combined.txt.gz")

        # Declare variables.
        self.eqtl_probes = None

    def start(self):
        print("Starting combining eQTL probe files.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.eqtl_probes = load_dataframe(inpath=self.outpath, header=0,
                                              index_col=False)
        else:
            # Load each GTE file.
            print("Loading eQTLprobes files.")
            self.eqtl_probes = self.combine_files()
            self.save()

    def combine_files(self):
        combined = None
        for i in range(1, self.n_iterations+1):
            infile = os.path.join(self.indir, self.iter_dirname + str(i),
                                  self.in_filename)
            df = load_dataframe(inpath=infile, header=0, index_col=False)
            df["Iteration"] = i
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)

        # Remove duplicate entries.
        combined.drop_duplicates(inplace=True)

        return combined

    def save(self):
        save_dataframe(df=self.eqtl_probes, outpath=self.outpath,
                       index=False, header=True)

    def clear_variables(self):
        self.indir = None
        self.iter_dirname = None
        self.n_iterations = None
        self.in_filename = None
        self.force = None

    def get_outpath(self):
        return self.outpath

    def get_eqtlprobes(self):
        return self.eqtl_probes

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Iteration directory: {}".format(self.iter_dirname))
        print("  > N. Iterations: {}".format(self.n_iterations))
        print("  > Input filename: {}".format(self.in_filename))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
