#!/usr/bin/env python3

"""
File:         optimize_covariates.py
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
from __future__ import print_function
from pathlib import Path
import itertools
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Optimize Covariates"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
./optimize_covariates.py -g /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/create_matrices/genotype_table.txt -e /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/create_matrices/expression_table.txt.gz -c /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/perform_deconvolution/deconvolution_table.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.genotype_path = getattr(arguments, 'genotype')
        self.expression_path = getattr(arguments, 'expression')
        self.covariates_path = getattr(arguments, 'covariates')

        # Set variables.
        self.outdir = str(Path(__file__).parent.parent)

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
        parser.add_argument("-g",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-e",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-c",
                            "--covariates",
                            type=str,
                            required=True,
                            help="The path to the covariates matrix")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        # Load the dataframes.
        print("Loading dataframes.")
        geno_df = self.load_dataframe(self.genotype_path)
        expr_df = self.load_dataframe(self.expression_path)
        cov_df = self.load_dataframe(self.covariates_path, nrows=None)

        # Validating input.
        print("Validating input.")
        self.validate_dataframes(geno_df, expr_df, cov_df)

        # Correct gene expression for cohorts.


        print(geno_df)
        print(expr_df)
        print(cov_df)

    @staticmethod
    def load_dataframe(inpath, header=0, index_col=0, sep="\t", low_memory=True,
                       nrows=100, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} with shape: {}".format(os.path.basename(inpath),
                                                             df.shape))
        return df

    @staticmethod
    def validate_dataframes(geno_df, expr_df, cov_df):
        dfs = [geno_df, expr_df, cov_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.columns.identical(b.columns):
                print("Order of samples are not identical.")
                exit()

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype datafile: {}".format(self.genotype_path))
        print("  > Expression datafile: {}".format(self.expression_path))
        print("  > Covariates datafile: {}".format(self.covariates_path))


if __name__ == '__main__':
    m = main()
    m.start()
