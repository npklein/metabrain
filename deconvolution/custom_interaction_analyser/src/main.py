"""
File:         main.py
Created:      2020/03/23
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
import os

# Third party imports.
import numpy as np
import pandas as pd
import pandas.io.common
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import scipy.stats

# Local application imports.
from general.local_settings import LocalSettings
from general.utilities import prepare_output_dir, p_value_to_symbol


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, force, verbose):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param force: boolean, whether or not to force to redo each step.
        :param verbose: boolean, whether or not to print each step.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.eqtl_inpath = settings.get_setting("eqtl_datafile")
        self.geno_inpath = settings.get_setting("genotye_datafile")
        self.expr_inpath = settings.get_setting("expression_datafile")
        self.cov_inpath = settings.get_setting("covariates_datafile")
        self.tech_covs = settings.get_setting("technical_covariates")
        self.force = force
        self.verbose = verbose

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Global variables.
        self.chunk_size = 10

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        print("Loading covariates")
        tech_cov, cov_df = self.get_covariates()
        n_samples = cov_df.shape[1]
        n_covariates = cov_df.shape[0]
        n_tech_covariates = tech_cov.shape[0]
        print("\tNumber of samples: {}".format(n_samples))
        print("\tNumber of covariates: {}".format(n_covariates))
        print("\tNumber of technical covariates: {}\n".format(n_tech_covariates))

        print("Processing")
        index = 0
        while True:
            try:
                eqtl_df, geno_df, expr_df = self.get_data(index)
            except pandas.io.common.EmptyDataError:
                break

            # Replace -1 with NaN.
            geno_df.replace(-1, np.nan, inplace=True)

            # loop over eQTLs.
            for i in range(self.chunk_size):
                # Get the missing genotype indices.
                indices = np.arange(n_samples)
                sample_mask = indices[~geno_df.iloc[i, :].isnull().values]
                self.print_string("Sample mask:\tLength:{}\tMin:{}"
                                  "\tMax:{}".format(len(sample_mask),
                                                    min(sample_mask),
                                                    max(sample_mask)))

                # Subset.
                genotype = geno_df.iloc[i, sample_mask].copy()
                expression = expr_df.iloc[i, sample_mask].copy()
                technical_covs = tech_cov.iloc[:, sample_mask].copy()
                covariates = cov_df.iloc[:, sample_mask].copy()

                # Check if SNP index is identical.
                if (genotype.index.name != expression.index.name):
                    print("Indices do not match.")
                    exit()

                # Create the null model.
                null_matrix = technical_covs.copy().mul(genotype, axis=1)
                y_null = self.multivariate_linear_regression(null_matrix.T,
                                                             expression)

                #  Calculate the statistics of the model.
                mse_null = mean_squared_error(y_null, expression)
                r2_null = r2_score(y_null, expression)
                self.print_string("Null model:\tMSE: {:.2f}\t"
                                  "R2: {:.2f}".format(mse_null, r2_null))

                # Loop over the covariates.
                for j in range(covariates.shape[0]):
                    # Append to covariate.
                    name = covariates.index[j]
                    print("Analyzing covariate: {}".format(name))
                    cov_of_interest = covariates.iloc[j, :] * genotype
                    covariates_matrix = null_matrix.copy()
                    covariates_matrix.loc[name, :] = cov_of_interest

                    # Create the alternative model.
                    alt_matrix = covariates_matrix.mul(genotype, axis=1)
                    y_alt = self.multivariate_linear_regression(alt_matrix.T,
                                                                expression)

                    #  Calculate the statistics of the model.
                    mse_alt = mean_squared_error(y_alt, expression)
                    r2_alt = r2_score(y_alt, expression)
                    self.print_string("  Alt model:\tMSE: {:.2f}\t"
                                      "R2: {:.2f}".format(mse_alt, r2_alt))

                    # Calculate the difference.
                    r2_diff = r2_alt - r2_null
                    p_value = scipy.stats.f.cdf(r2_diff, dfn=1, dfd=1)
                    self.print_string(
                        "  Model Comparison:\tDelta R2: {:.2f}\t"
                        "p-value: {:.2e} [{}]".format(r2_diff, p_value,
                                                      p_value_to_symbol(p_value)))

            # Increment index and check if end is reached.
            if len(eqtl_df.index) != self.chunk_size:
                break
            index += 1

            exit()

    def get_covariates(self):
        cov_df = pd.read_csv(self.cov_inpath,
                             sep="\t",
                             header=0,
                             index_col=0)

        tech_cov = cov_df.loc[self.tech_covs, :]
        cov_df.drop(self.tech_covs, inplace=True)

        return tech_cov, cov_df

    def get_data(self, index):
        self.print_string("\tLoading index {}-{}".format(index * self.chunk_size,
                                                         (index+1) * self.chunk_size))
        eqtl_df = pd.read_csv(self.eqtl_inpath,
                              sep="\t",
                              header=0,
                              skiprows=index * self.chunk_size,
                              nrows=self.chunk_size)

        geno_df = pd.read_csv(self.geno_inpath,
                              sep="\t",
                              header=0,
                              index_col=0,
                              skiprows=index * self.chunk_size,
                              nrows=self.chunk_size)

        expr_df = pd.read_csv(self.expr_inpath,
                              sep="\t",
                              header=0,
                              index_col=0,
                              skiprows=index * self.chunk_size,
                              nrows=self.chunk_size)

        return eqtl_df, geno_df, expr_df

    def multivariate_linear_regression(self, X, y):
        self.print_string("X = {}, y = {}".format(X.shape, y.shape))
        regressor = LinearRegression()
        regressor.fit(X, y)
        y_pred = regressor.predict(X)
        return y_pred

    def print_string(self, string):
        if self.verbose:
            print(string)

    def print_arguments(self):
        print("Arguments:")
        print("  > EQTL datafile: {}".format(self.eqtl_inpath))
        print("  > Genotype datafile: {}".format(self.geno_inpath))
        print("  > Expression datafile: {}".format(self.expr_inpath))
        print("  > Covariates datafile: {}".format(self.cov_inpath))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > Force: {}".format(self.force))
        print("  > Verbose: {}".format(self.verbose))
        print("  > Chunk size: {}".format(self.chunk_size))
        print("")
