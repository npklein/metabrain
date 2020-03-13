#!/usr/bin/env python3

"""
File:         combine_eQTLprobes.py
Created:      2020/02/25
Last Changed: 2020/03/03
Author:       M.Vochteloo

Copyright (C) 2019 M.Vochteloo

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
import glob
import os

# Third party imports.
import pandas as pd

# Local application imports.


# Metadata.
__program__ = "Combine eQTLprobes Files"
__author__ = "M. Vochteloo"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class Main:
    """
    Main class of the program.
    """
    def __init__(self, indir, in_filename, n_iter):
        """
        The initializer for the main class.

        :param indir: string, the input directory containing GTE files.
        :param in_filename: string, the input filename.
        :param n_iter: int, the number of iterations to include.
        """
        self.indir = indir
        self.in_filename = in_filename
        self.n_iter = n_iter

        outdir = os.path.join(os.getcwd(), 'output')
        self.outpath = os.path.join(outdir, in_filename + "_combined.txt.gz")

        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def start(self):
        """
        Main method for the main class. Does all the work.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Input filename: {}".format(self.in_filename))
        print("  > N iterations: {}".format(self.n_iter))
        print("  > Output path: {}".format(self.outpath))
        print("")

        # Combine each infile into a dataframe.
        print("Loading infiles.")
        combined = None
        for i in range(1, self.n_iter+1):
            filename = os.path.join(self.indir, "Iteration" + str(i),
                                    self.in_filename + ".txt.gz")
            df = pd.read_csv(filename, sep="\t", header=0)
            print("\t[{}] File: {}, shape: {}".format(i, filename, df.shape))
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)
        print("")
        print("\tShape: {}".format(combined.shape))

        # Remove duplicate entries.
        print("Filtering on duplicates.")
        combined.drop_duplicates(inplace=True)
        print("\tShape: {}".format(combined.shape))

        # write outfile.
        combined.to_csv(self.outpath, sep="\t", index=False, compression='gzip')


if __name__ == "__main__":
    # Define main variables.
    INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output",
                         "2019-11-06-FreezeTwoDotOne",
                         "2020-02-18-eqtls",
                         "cortex-cis-EURandAFR-iterative",
                         )

    IN_FILENAME = "eQTLProbesFDR0.05-ProbeLevel"
    N_ITER = 4

    # Start the program.
    MAIN = Main(indir=INDIR,
                in_filename=IN_FILENAME,
                n_iter=N_ITER)
    MAIN.start()
