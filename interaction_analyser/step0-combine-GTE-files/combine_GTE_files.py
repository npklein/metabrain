#!/usr/bin/env python3

"""
File:         combine_GTE_files.py
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
__program__ = "Combine GTE files"
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
    def __init__(self, indir, out_filename):
        """
        The initializer for the main class.

        :param indir: string, the input directory containing GTE files.
        :param out_filename: string, the output filename.
        """
        self.indir = indir
        outdir = os.path.join(os.getcwd(), 'output')
        self.outpath = os.path.join(outdir, out_filename)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def start(self):
        """
        Main method for the main class. Does all the work.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output path: {}".format(self.outpath))
        print("")

        # Combine each GTE files into a dataframe.
        print("Loading GTE files.")
        combined = None
        for i, infile in enumerate(glob.glob(os.path.join(self.indir, "GTE-*.txt"))):
            df = pd.read_csv(infile, sep="\t", header=None)
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)
            print("\t[{}] File: {}, shape: {}".format(i, infile, df.shape))
        print("")
        print("\tShape: {}".format(combined.shape))

        # Remove duplicate entries.
        print("Filtering on duplicates.")
        combined.drop_duplicates(inplace=True)
        print("\tShape: {}".format(combined.shape))

        # write outfile.
        combined.to_csv(self.outpath, sep="\t", index=False, header=False,
                        compression='gzip')


if __name__ == "__main__":
    # Define main variables.
    INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output",
                         "2019-11-06-FreezeTwoDotOne",
                         "2020-02-18-eqtls",
                         "cortex-cis-EURandAFR-iterative",
                         )

    OUT_FILENAME = "GTE_combined.txt.gz"

    # Start the program.
    MAIN = Main(indir=INDIR, out_filename=OUT_FILENAME)
    MAIN.start()
