#!/usr/bin/env python3

"""
File:         add_ct_mediated_info_to_matrix.py
Created:      2020/10/28
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
__program__ = "Add Cell Type Mediated Info to Matrix"
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
        self.matrix_path = getattr(arguments, 'matrix')
        self.skip_rows = getattr(arguments, 'skip_rows')
        self.probe_id = getattr(arguments, 'probe_id')
        self.snp_id = getattr(arguments, 'snp_id')
        self.decon_path = getattr(arguments, 'decon')

        # Set variables.
        self.outfile = os.path.join(Path(self.matrix_path).parent, self.get_basename(self.matrix_path) + "_withCTMediationScores.txt.gz")

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
        parser.add_argument("-m",
                            "--matrix",
                            type=str,
                            required=True,
                            help="The path to the matrix.")
        parser.add_argument("-s",
                            "--skip_rows",
                            type=int,
                            default=0,
                            help="The number of rows to skip. Default: 0")
        parser.add_argument("-probe",
                            "--probe_id",
                            type=str,
                            required=True,
                            help="The probe column name.")
        parser.add_argument("-snp",
                            "--snp_id",
                            nargs="*",
                            type=str,
                            required=True,
                            help="The snp column name(s).")
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
            if part not in ["txt", "tsv", "gz", "xlsx"]:
                new_basename.append(part)

        return ".".join(new_basename)

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data.")
        matrix_df = self.load_file(self.matrix_path, skiprows=self.skip_rows,
                                   index_col=None)
        print(matrix_df)
        decon_df = self.load_file(self.decon_path)

        print("### Step2 ###")
        print("Calculate BH-FDR.")
        decon_fdr_df = self.bh_correct(decon_df)
        print(decon_fdr_df)

        print("### Step3 ###")
        print("Preprocessing.")

        snp_id_col = "SNPName"
        if len(self.snp_id) > 1:
            matrix_snp_names = []
            for _, row in matrix_df.iterrows():
                snp_name = []
                for col in self.snp_id:
                    snp_name.append(str(row[col]))
                matrix_snp_names.append(":".join(snp_name))
            matrix_df[snp_id_col] = matrix_snp_names
        else:
            snp_id_col = self.snp_id[0]
        print(matrix_df[[snp_id_col]])

        probe_names = []
        snp_names = []
        for index in decon_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append(index.split("_")[1])
        decon_fdr_df["ProbeName"] = probe_names
        decon_fdr_df["SNPName"] = snp_names
        print(decon_fdr_df)

        print(matrix_df.loc[matrix_df[self.probe_id] == "ENSG00000019186.10", :])
        print(decon_fdr_df.loc[decon_fdr_df["ProbeName"] == "ENSG00000019186.10", :])

        print("### Step4 ###")
        print("Merge.")
        merged_df = matrix_df.merge(decon_fdr_df,
                                    left_on=[self.probe_id, snp_id_col],
                                    right_on=["ProbeName", "SNPName"])
        merged_df = merged_df.drop(['ProbeName', 'SNPName'], axis=1)

        self.save_file(merged_df,
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
        print("  > Matrix path: {}".format(self.matrix_path))
        print("  > Skip rows: {}".format(self.skip_rows))
        print("  > Probe ID: {}".format(self.probe_id))
        print("  > SNP ID: {}".format(self.snp_id))
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > Output file: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
