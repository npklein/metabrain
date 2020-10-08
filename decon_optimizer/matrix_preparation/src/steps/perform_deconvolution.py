"""
File:         perform_deconvolution.py
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
import numpy as np
import pandas as pd
from scipy.optimize import nnls
import statsmodels.api as sm

# Local application imports.
from matrix_preparation.src.utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class PerformDeconvolution:
    def __init__(self, settings, sign_file, sign_df, sign_expr_file,
                 sign_expr_df, cohort_df, force, outdir):
        self.sign_file = sign_file
        self.sign_df = sign_df
        self.sign_expr_file = sign_expr_file
        self.sign_expr_df = sign_expr_df
        self.cohort_df = cohort_df
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'perform_deconvolution')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.outpath = os.path.join(self.outdir, "deconvolution_table.txt.gz")

        # Create empty variable.
        self.decon_df = None

    def start(self):
        print("Starting deconvolution.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.decon_df = load_dataframe(inpath=self.outpath,
                                                header=0, index_col=0)
        else:
            self.decon_df = self.perform_deconvolution()
            self.save()

    def perform_deconvolution(self):
        if self.sign_df is None:
            # Load the celltype profile file.
            print("Loading cell type profile matrix.")
            self.sign_df = load_dataframe(self.sign_file,
                                          header=0,
                                          index_col=0)

        if self.sign_expr_df is None:
            # Load the celltype expression file.
            print("Loading cell type expression matrix.")
            self.sign_expr_df = load_dataframe(self.sign_expr_file,
                                               header=0,
                                               index_col=0)

        # Filter uninformative genes from the signature matrix.
        sign_df = self.filter(self.sign_df)

        # Subset and reorder.
        sign_df, expr_df, cohort_df = self.subset(sign_df,
                                                  self.sign_expr_df,
                                                  self.cohort_df)

        # Transform.
        sign_df = self.perform_log2_transform(sign_df)

        # Shift the data to be positive.
        print("Shifting data to be positive")
        if sign_df.values.min() < 0:
            print("\tSignature matrix is shifted.")
            sign_df = self.perform_shift(sign_df)

        if expr_df.values.min() < 0:
            print("\tExpression matrix is shifted.")
            expr_df = self.perform_shift(expr_df)

        print("Signature shape: {}".format(sign_df.shape))
        print("Expression shape: {}".format(expr_df.shape))

        # Perform deconvolution per sample.
        print("Performing partial deconvolution.")
        decon_data = []
        residuals_data = []
        for _, sample in expr_df.T.iterrows():
            proportions, rnorm = self.nnls(sign_df, sample)
            decon_data.append(proportions)
            residuals_data.append(rnorm)

        decon_df = pd.DataFrame(decon_data,
                                index=expr_df.columns,
                                columns=["{}NNLS_{}".format(*x.split("_")) for x in sign_df.columns])
        residuals_df = pd.Series(residuals_data, index=expr_df.columns)

        print("Estimated weights:")
        print(decon_df)
        print(decon_df.mean(axis=0))

        # Make the weights sum up to 1.
        decon_df = self.sum_to_one(decon_df)
        print("Estimated proportions:")
        print(decon_df)
        print(decon_df.mean(axis=0))

        # Calculate the average residuals.
        print(residuals_df)
        print("Average residual: {:.2f}".format(residuals_df.mean()))

        return decon_df

    @staticmethod
    def filter(df, cutoff=0):
        tmp = df.copy()
        return tmp.loc[(tmp.std(axis=1) != 0) & (tmp.max(axis=1) > cutoff), :]

    @staticmethod
    def subset(raw_signature, raw_expression, raw_cohorts):
        gene_overlap = np.intersect1d(raw_signature.index, raw_expression.index)
        sign_df = raw_signature.loc[gene_overlap, :]
        expr_df = raw_expression.loc[gene_overlap, :]

        sign_df.index.name = "-"
        expr_df.index.name = "-"

        if not sign_df.index.equals(expr_df.index):
            print("Invalid gene order")
            exit()

        cohorts_df = None
        if raw_cohorts is not None:
            sample_overlap = np.intersect1d(expr_df.columns, raw_cohorts.index)
            expr_df = expr_df.loc[:, sample_overlap]
            cohorts_df = raw_cohorts.loc[sample_overlap, :]

            if not expr_df.columns.equals(cohorts_df.index):
                print("Invalid sample order")
                exit()

        return sign_df, expr_df, cohorts_df

    @staticmethod
    def perform_log2_transform(df):
        tmp = df.copy()
        tmp = tmp + abs(tmp.values.min()) + 1
        return np.log2(tmp)

    @staticmethod
    def perform_shift(df):
        return df + abs(df.values.min())

    @staticmethod
    def nnls(A, b):
        return nnls(A, b)

    @staticmethod
    def sum_to_one(X):
        return X.divide(X.sum(axis=1), axis=0)

    def save(self):
        save_dataframe(df=self.decon_df, outpath=self.outpath,
                       index=True, header=True)

    def clear_variables(self):
        self.sign_file = None
        self.sign_df = None
        self.sign_expr_file = None
        self.sign_expr_df = None
        self.cohort_df = None
        self.force = None

    def get_decon_df(self):
        return self.decon_df

    def print_arguments(self):
        print("Arguments:")
        if self.sign_df is not None:
            print("  > Signature: {}".format(self.sign_df.shape))
        else:
            print("  > Signature input file: {}".format(self.sign_file))
        if self.sign_expr_df is not None:
            print("  > Signature expression: {}".format(self.sign_expr_df.shape))
        else:
            print("  > Signature input path: {}".format(self.sign_expr_file))
        print("  > Deconvolution output file: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
