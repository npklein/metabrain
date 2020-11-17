#!/usr/bin/env python3

"""
File:         create_decon_fdr_df.py
Created:      2020/11/16
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
from __future__ import print_function
from pathlib import Path
import argparse
import os

# Third party imports.
import pandas as pd
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Create Decon FDR DF"
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
        self.decon_path = getattr(arguments, 'decon')
        self.outfile = getattr(arguments, 'outfile')

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'output')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # Define the output file.
        self.outpath = os.path.join(outdir, self.outfile + ".txt.gz")

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            default="deconvolutionResults_FDR",
                            help="The name of the ouput file. Default: "
                                 "'deconvolutionResults_FDR'")

        return parser.parse_args()

    def start(self):
        print("### Step1 ###")
        print("Loading data.")
        decon_df = pd.read_csv(self.decon_path,
                               sep="\t",
                               header=0,
                               index_col=0)
        print("\tDeconvolution results data frame: {}".format(decon_df.shape))

        print("### Step2 ###")
        print("Calculate BH-FDR.")
        decon_fdr_df = self.bh_correct(decon_df)
        cols = decon_fdr_df.columns.tolist()

        print("Preprocessing deconvolution data frame")
        probe_names = []
        snp_names = []
        for index in decon_fdr_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        decon_fdr_df["ProbeName"] = probe_names
        decon_fdr_df["SNPName"] = snp_names

        print("Reordering columns")
        cols.sort()
        cols = ["SNPName", "ProbeName"] + cols
        decon_fdr_df = decon_fdr_df.loc[:, cols]
        print(decon_fdr_df)

        print("Saving data frame.")
        decon_fdr_df.to_csv(self.outpath, sep="\t", compression="gzip",
                            header=True, index=False)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(self.outpath),
                                      decon_fdr_df.shape))

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        fdr_df = pd.DataFrame(data, index=indices, columns=df.index)

        return fdr_df.T

    def print_arguments(self):
        print("Arguments:")
        print("  > Decon file: {}".format(self.decon_path))
        print("  > Output file: {}".format(self.outfile))
        print("  > Outpath {}".format(self.outpath))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
