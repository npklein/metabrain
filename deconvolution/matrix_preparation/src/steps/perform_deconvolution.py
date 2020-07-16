"""
File:         perform_deconvolution.py
Created:      2020/04/08
Last Changed: 2020/07/15
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

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe, save_dataframe


class PerformDeconvolution:
    def __init__(self, settings, profile_file, profile_df, ct_expr_file,
                 ct_expr_df, force, outdir):
        """
        The initializer for the class.


        :param settings: string, the settings.
        :param profile_file: string, the datafile contaioning the celltype
                             profile.
        :param profile_df: DataFrame, the celltype profile.
        :param ct_expr_file: string, the datafile containing expression
                             of the celltype profiles.
        :param ct_expr_df: string, the celltype expression.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.profile_file = profile_file
        self.profile_df = profile_df
        self.ct_expr_file = ct_expr_file
        self.ct_expr_df = ct_expr_df
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'perform_deconvolution')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.outpath = os.path.join(self.outdir, "deconvolution_table.txt.gz")

        # Create empty variable.
        self.deconvolution = None

    def start(self):
        print("Starting deconvolution.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.deconvolution = load_dataframe(inpath=self.outpath,
                                                header=0, index_col=0)
        else:
            self.deconvolution = self.perform_deconvolution()
            self.save()

    def perform_deconvolution(self):
        if self.profile_df is None:
            # Load the celltype profile file.
            print("Loading cell type profile matrix.")
            self.profile_df = load_dataframe(self.profile_file,
                                             header=0, index_col=0)

        if self.ct_expr_df is None:
            # Load the celltype expression file.
            print("Loading cell type expression matrix.")
            self.ct_expr_df = load_dataframe(self.ct_expr_file,
                                             header=0, index_col=0)

        # Filter uninformative genes from the signature matrix.
        prof_df = self.filter(self.profile_df)

        # Subset and reorder.
        prof_df, expr_df = self.subset(prof_df, self.ct_expr_df)

        # Transform.
        prof_df = self.perform_log2_transform(prof_df)

        # Shift the data to be positive.
        print("Shifting data to be positive")
        if prof_df.values.min() < 0:
            prof_df = self.perform_shift(prof_df)

        if expr_df.values.min() < 0:
            expr_df = self.perform_shift(expr_df)

        print("Profile shape: {}".format(prof_df.shape))
        print("Expression shape: {}".format(expr_df.shape))

        # Perform deconvolution per sample.
        print("Performing partial deconvolution.")
        decon_data = []
        residuals_data = []
        for _, sample in expr_df.T.iterrows():
            proportions, rnorm = self.nnls(prof_df, sample)
            decon_data.append(proportions)
            residuals_data.append(rnorm)

        decon_df = pd.DataFrame(decon_data,
                                index=expr_df.columns,
                                columns=["{}NNLS_{}".format(*x.split("_")) for x in prof_df.columns])
        residuals_df = pd.Series(residuals_data, index=expr_df.columns)

        print("Estimated weights:")
        print(decon_df)
        print(decon_df.mean(axis=0))

        # Make the weights sum up to 1.
        decon_df = self.sum_to_one(decon_df)
        print("Estimated proportions:")
        print(decon_df)

        # Calculate the average residuals.
        print(residuals_df)
        print("Average residual: {:.2f}".format(residuals_df.mean()))

        return decon_df

    @staticmethod
    def filter(df, cutoff=0):
        tmp = df.copy()
        return tmp.loc[(tmp.std(axis=1) != 0) & (tmp.max(axis=1) > cutoff), :]

    @staticmethod
    def subset(raw_signature, raw_expression):
        overlap = np.intersect1d(raw_signature.index, raw_expression.index)
        sign_df = raw_signature.loc[overlap, :]
        expr_df = raw_expression.loc[overlap, :]

        sign_df.index.name = "-"
        expr_df.index.name = "-"

        if not sign_df.index.equals(expr_df.index):
            print("Invalid order")
            exit()

        return sign_df, expr_df

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
        save_dataframe(df=self.deconvolution, outpath=self.outpath,
                       index=True, header=True)

    def clear_variables(self):
        self.profile_file = None
        self.profile_df = None
        self.ct_expr_file = None
        self.ct_expr_df = None
        self.force = None

    def get_deconvolution(self):
        return self.deconvolution

    def print_arguments(self):
        print("Arguments:")
        if self.profile_df is not None:
            print("  > Celltype profile: {}".format(self.profile_df.shape))
        else:
            print("  > Celltype profile input file: {}".format(self.profile_file))
        if self.ct_expr_df is not None:
            print("  > Celltype expression: {}".format(self.ct_expr_df.shape))
        else:
            print("  > Celltype expression input path: {}".format(self.ct_expr_file))
        print("  > Deconvolution output file: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
