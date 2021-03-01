"""
File:         perform_deconvolution.py
Created:      2020/04/08
Last Changed: 2020/09/15
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
        self.sample_cohort_file = settings["sample_cohort_datafile"]
        self.sample_id = settings["sample_cohort_identifiers"]["sample"]
        self.cohort_id = settings["sample_cohort_identifiers"]["cohort"]
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

        print("Loading sample cohort matrix.")
        sample_cohort_df = load_dataframe(self.sample_cohort_file,
                                          header=0,
                                          index_col=None)

        # Correct for cohort effects.
        cohort_df = self.create_cohort_df(list(self.ct_expr_df.columns),
                                          sample_cohort_df,
                                          self.sample_id,
                                          self.cohort_id)

        # Filter uninformative genes from the signature matrix.
        prof_df = self.filter(self.profile_df)

        # Subset and reorder.
        prof_df, expr_df, cohort_df = self.subset(prof_df,
                                                  self.ct_expr_df,
                                                  cohort_df)

        # Correct for cohorts.
        expr_df = self.cohort_correction(expr_df, cohort_df)

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
    def create_cohort_df(samples, sample_cohort_df, sample_id, cohort_id):
        sample_to_cohort_dict = dict(zip(sample_cohort_df.loc[:, sample_id], sample_cohort_df.loc[:, cohort_id]))

        cohort_to_sample_dict = {}
        for sample in samples:
            cohort = sample_to_cohort_dict[sample]
            if cohort in cohort_to_sample_dict.keys():
                value = cohort_to_sample_dict[cohort]
                value.append(sample)
                cohort_to_sample_dict[cohort] = value
            else:
                cohort_to_sample_dict[cohort] = [sample]

        cohort_df = pd.DataFrame(0,
                                 index=sample_to_cohort_dict.keys(),
                                 columns=cohort_to_sample_dict.keys())
        for cohort in cohort_df.columns:
            cohort_df.loc[cohort_df.index.isin(cohort_to_sample_dict[cohort]), cohort] = 1

        return cohort_df

    @staticmethod
    def cohort_correction(raw_expression, cohorts):
        new_expression_data = []
        for i, (index, expression) in enumerate(raw_expression.iterrows()):
            if (i % 100 == 0) or (i == (raw_expression.shape[0] - 1)):
                print("\tProcessing {}/{} "
                      "[{:.2f}%]".format(i,
                                         (raw_expression.shape[0] - 1),
                                         (100 / (raw_expression.shape[0] - 1)) * i))

            if expression.index.equals(cohorts.index):
                ols = sm.OLS(expression, cohorts)
                try:
                    ols_result = ols.fit()
                    residuals = ols_result.resid.values
                    corrected_expression = expression.mean() + residuals
                    new_expression_data.append(corrected_expression.tolist())

                except np.linalg.LinAlgError as e:
                    print("\t\tError: {}".format(e))

        new_expression_df = pd.DataFrame(new_expression_data,
                                         index=raw_expression.index,
                                         columns=raw_expression.columns)

        return new_expression_df

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
