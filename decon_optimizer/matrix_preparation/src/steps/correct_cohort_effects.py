"""
File:         correct_cohort_effects.py
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
from matrix_preparation.src.utilities import prepare_output_dir, check_file_exists, save_dataframe


class CorrectCohortEffects:
    def __init__(self, settings, cohort_df, expr_df, sign_expr_df, force,
                 outdir):
        self.cohort_df = cohort_df
        self.expr_df = expr_df
        self.sign_expr_df = sign_expr_df
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'correct_cohort_effects')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.expr_cc_outpath = os.path.join(self.outdir, "expression_table_cc.txt.gz")
        self.sign_expr_cc_outpath = os.path.join(self.outdir, "signature_expr_table_cc.txt.gz")

        # Create empty variable.
        self.expr_cc_df = None
        self.sign_expr_cc_df = None

    def start(self):
        print("Starting deconvolution.")
        self.print_arguments()

        if not check_file_exists(self.expr_cc_outpath) or self.force:
            self.expr_cc_df = self.cohort_correction(self.expr_df, self.cohort_df)
            save_dataframe(df=self.expr_cc_df, outpath=self.expr_cc_outpath,
                           index=True, header=True)

        if not check_file_exists(self.sign_expr_cc_outpath) or self.force:
            self.sign_expr_df = self.cohort_correction(self.sign_expr_df, self.cohort_df)
            save_dataframe(df=self.expr_cc_df, outpath=self.sign_expr_cc_outpath,
                           index=True, header=True)

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

    def clear_variables(self):
        self.cohort_df = None
        self.expr_df = None
        self.sign_expr_df = None
        self.force = None

    def get_expr_cc_df(self):
        return self.expr_cc_df

    def get_sign_expr_cc_df(self):
        return self.sign_expr_cc_df

    def print_arguments(self):
        print("Arguments:")
        print("  > Cohort input shape: {}".format(self.cohort_df.shape))
        print("  > Expression input shape: {}".format(self.expr_df.shape))
        print("  > Signature expression input shape: {}".format(self.sign_expr_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Force: {}".format(self.force))
