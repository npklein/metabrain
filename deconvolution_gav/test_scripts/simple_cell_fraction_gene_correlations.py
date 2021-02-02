#!/usr/bin/env python3

"""
File:         simple_cell_fraction_gene_correlations.py
Created:      2021/02/02
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
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Simple Cell Fractions Gene Correlations"
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

"""
Syntax:
./simple_cell_fraction_gene_correlations.py -cf ../matrix_preparation/cortex_eur_trans/perform_deconvolution/deconvolution_table.txt -ge ../matrix_preparation/cortex_eur_trans/create_matrices/expression_table.txt.gz -op cortex_eur_trans_
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.ge_path = getattr(arguments, 'gene_expression')
        self.outfile_prefix = getattr(arguments, 'outfile_prefix')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "simple_gene_cellfraction_corr")

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
                            help="show program's version number and exit")
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fracttions matrix")
        parser.add_argument("-ge",
                            "--gene_expression",
                            type=str,
                            required=True,
                            help="The path to the gene expression matrix")
        parser.add_argument("-op",
                            "--outfile_prefix",
                            type=str,
                            required=False,
                            default="",
                            help="A prefix for the output files.")

        return parser.parse_args()

    def start(self):
        print("Loading cell fractions matrix.")
        cf_df = self.load_file(self.cf_path, low_memory=False)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cf_path),
                                      cf_df.shape))

        print("Loading gene expression matrix.")
        ge_df = self.load_file(self.ge_path)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.ge_path),
                                      ge_df.shape))

        print("Pre-process")
        if not cf_df.index.equals(ge_df.columns):
            print("The expressiom file columns do not match the cell "
                  "type fractions file.")
            exit()
        cf_df = cf_df.T
        ge_df = ge_df.drop_duplicates()

        print(cf_df)
        print(ge_df)

        print("Correlate data.")
        corr_coefficients = []
        pvalues = []
        for cell_type, cf_data in cf_df.iterrows():
            print("\tProcessing '{}'".format(cell_type))
            ct_coefficients = []
            ct_pvalues = []
            for i, (gene, expr_data) in enumerate(ge_df.iterrows()):
                coef, p = stats.spearmanr(cf_data, expr_data)
                ct_coefficients.append(coef)
                ct_pvalues.append(p)
            corr_coefficients.append(ct_coefficients)
            pvalues.append(ct_pvalues)

        corr_coef_df = pd.DataFrame(corr_coefficients,
                                    index=cf_df.index,
                                    columns=ge_df.index)
        print(corr_coef_df)
        pvalue_df = pd.DataFrame(pvalues,
                                 index=cf_df.index,
                                 columns=ge_df.index)
        print(pvalue_df)

        # Save.
        corr_coef_df.to_csv(os.path.join(self.outdir, '{}coefficients.txt.gz'.format(self.outfile_prefix)),
                            sep="\t", index=True, header=True, compression="gzip")
        pvalue_df.to_csv(os.path.join(self.outdir, '{}pvalues.txt.gz'.format(self.outfile_prefix)),
                         sep="\t", index=True, header=True, compression="gzip")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None, low_memory=True):
        if path.endswith(".pkl"):
            df = pd.read_pickle(path)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df


if __name__ == '__main__':
    m = main()
    m.start()
