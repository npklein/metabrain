#!/usr/bin/env python3

"""
File:         decon_qtl_parser.py
Created:      2020/11/20
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
__program__ = "Decon QTL Parser"
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
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.decon_path = getattr(arguments, 'decon')

        # Set variables.
        outdir = Path(__file__).parent.absolute()
        self.outfile = os.path.join(outdir, self.get_basename(self.decon_path) + "alleles_FDR_betas.txt.gz")

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-e",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix using for decon.")
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix.")

        return parser.parse_args()

    @staticmethod
    def get_basename(filepath):
        basename = os.path.basename(filepath)
        splitted_name = basename.split(".")
        new_basename = []
        for part in splitted_name:
            if part not in ["txt", "csv", "tsv", "gz", "xlsx"]:
                new_basename.append(part)

        return ".".join(new_basename)

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data.")
        eqtl_df = self.load_file(self.eqtl_path, index_col=None)
        print(eqtl_df)
        decon_df = self.load_file(self.decon_path)
        print(decon_df)

        print("### Step2 ###")
        print("Preprocessing.")
        probe_names = []
        snp_names = []
        for index in decon_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        decon_df["ProbeName"] = probe_names
        decon_df["SNPName"] = snp_names
        print(decon_df)

        print("### Step4 ###")
        print("Merge.")
        merged_df = eqtl_df.merge(decon_df,
                                  left_on=["ProbeName", "SNPName"],
                                  right_on=["ProbeName", "SNPName"])
        print(merged_df)

        print("### Step5 ###")
        print("Adding BH-FDR.")
        fdr_columns = []
        for col in merged_df.columns:
            if col.endswith("_pvalue"):
                colname = col.replace("_pvalue", "_FDR")
                fdr_columns.append(colname)
                merged_df[colname] = multitest.multipletests(merged_df.loc[:, col], method='fdr_bh')[1]

        print("### Step5 ###")
        print("Filtering columns.")
        beta_colnames = []
        beta_columns = []
        for col in merged_df.columns:
            if col.endswith(":GT"):
                colname = "_".join([col.split("_")[1].replace(":GT", ""), "Beta"])
                beta_colnames.append(colname)
                beta_columns.append(col)
        filtered_df = merged_df.loc[:, ["SNPName", "ProbeName", "HGNCName", "SNPType", "AlleleAssessed"] + beta_columns + fdr_columns]
        print(filtered_df)
        filtered_df.columns = ["SNPName", "ProbeName", "HGNCName", "SNPType", "AlleleAssessed"] + beta_colnames + fdr_columns
        print(filtered_df)

        self.save_file(filtered_df,
                       self.outfile,
                       index=False)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", "_InteractionFDR"))
        fdr_df = pd.DataFrame(data, index=indices, columns=df.index)

        return fdr_df.T

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print(
            "\tSaved dataframe: {} "
            "with shape: {}".format(os.path.basename(outpath),
                                    df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Decon rows: {}".format(self.decon_path))
        print("  > Output file: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
