#!/usr/bin/env python3

"""
File:         decon_eqtl_with_permutation_fdr.py
Created:      2021/07/12
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
import itertools
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.special import betainc

# Local application imports.


"""
Syntax:
./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/expression_table/2020-07-16-MetaBrainDeconQtlGenes.TMM.SampSelect.ZeroVarRemov.covRemoved.expAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt.gz
"""

# Metadata
__program__ = "Decon-eQTL with Permutation FDR"
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
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.missing_geno = getattr(arguments, 'missing_genotype')
        self.alpha = getattr(arguments, 'alpha')
        self.n_permutations = getattr(arguments, 'permutations')
        outdir = getattr(arguments, 'outdir')
        self.print_interval = 30
        self.nrows = 1

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent.absolute()), "decon_eqtl_with_permutation_fdr", outdir)
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
                            help="show program's version number and exit.")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts "
                                 "matrix")
        parser.add_argument("-m",
                            "--missing_genotype",
                            type=int,
                            required=False,
                            default=-1,
                            help="The genotype value that equals a missing "
                                 "value. Default: -1. Note: has to be int.")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="The significance cut-off. Default: 0.05.")
        parser.add_argument("-p",
                            "--permutations",
                            type=int,
                            default=1000,
                            help="The number of permutations to run.")
        parser.add_argument("-o",
                            "--outdir",
                            type=str,
                            default="output",
                            help="The name of the ouput directory. "
                                 "Default: 'output'")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### STEP 1 ###")
        print("Loading data.")
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, nrows=self.nrows)
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, nrows=self.nrows)
        cc_df = self.load_file(self.cc_path, header=0, index_col=0)

        # Transpose if need be.
        if cc_df.shape[0] == geno_df.shape[1] or cc_df.shape[0] == expr_df.shape[1]:
            print("\t  Transposing covariate matrix.")
            cc_df = cc_df.T

        # Remove columns with all nan.
        expr_df.dropna(axis='columns', how='all', inplace=True)

        print("\tValidating input.")
        self.validate_data(geno_df=geno_df,
                           expr_df=expr_df,
                           cc_df=cc_df)
        print("")

        ########################################################################

        # print(geno_df)
        # print(expr_df)
        # print(cc_df)
        
        print("### STEP 2 ###")
        print("Pre-processing data.")
        # Check if the requirements are met.
        if expr_df.min(axis=1).min() < 0:
            print("Error: expression matrix contains negative values.")
            exit()
        if cc_df.min(axis=1).min() < 0:
            print("Error: cell count matrix contains negative values.")
            exit()

        # Convert to numpy.
        geno_m = geno_df.to_numpy(np.float64)
        expr_m = expr_df.to_numpy(np.float64)
        cc_m = cc_df.to_numpy(np.float64)
        
        # Save properties.
        eqtl_indices = expr_df.index + "_" + geno_df.index
        cell_types_indices = cc_df.index.to_numpy(dtype=object)
        del geno_df, expr_df, cc_df
        
        # print info
        n_samples = geno_m.shape[1]
        n_eqtls = geno_m.shape[0]
        n_covariates = cc_m.shape[0]
        print("\tN-samples: {}".format(n_samples))
        print("\tN-eQTLs: {}".format(n_eqtls))
        print("\tN-covariates: {}".format(n_covariates))
        print("")

        ########################################################################

        print("### STEP 3 ###")
        print("Analyzing.")
        output_m = np.empty((n_eqtls, n_covariates * 3), dtype=np.float64)
        config_alt_m = np.empty((n_eqtls, n_covariates), dtype=np.dtype('u1'))
        config_null_dict = {}
        for cell_type in cell_types_indices:
            empty_m = np.empty((n_eqtls, n_covariates))
            empty_m[:] = np.nan
            config_null_dict[cell_type] = empty_m
        alt_model_configs = list(itertools.product([True, False], repeat=n_covariates))
        null_model_configs = list(itertools.product([True, False], repeat=n_covariates - 1))
        start_time = int(time.time())
        last_print_time = None
        for row_index in range(n_eqtls):
            # Get the genotype.
            genotype = geno_m[row_index, :]

            # Construct the mask to remove missing values.
            mask = genotype != self.missing_geno
            n = np.sum(mask)

            # Select the values we use from genotype and expression.
            y = expr_m[row_index, mask]
            genotype = genotype[mask]

            # Model the alternative matrix (with the interaction term).
            config_alt, X_alt, betas_alt, rss_alt = \
                self.find_best_config(configs=alt_model_configs,
                                      n_samples=n,
                                      n_covariates=n_covariates,
                                      cc_m=cc_m[:, mask],
                                      genotype=genotype,
                                      expression=y
                                      )

            # Calculate the RSS just like Decon-eQTL does. This is not
            # how it was intended to work!
            # TODO delete this
            full_y = expr_m[row_index, :]
            alt_y_hat = np.zeros_like(full_y)
            alt_y_hat[mask] = self.predict(X=X_alt, betas=betas_alt)
            rss_alt = self.calc_rss(y=full_y, y_hat=alt_y_hat)

            # Save the alternative model stats.
            config_alt_m[row_index, :] = config_alt
            flip_array = np.vectorize({True: -1, False: 1}.get)(config_alt)
            output_m[row_index, n_covariates:2 * n_covariates] = betas_alt[:n_covariates]
            output_m[row_index, 2 * n_covariates:] = betas_alt[n_covariates:] * flip_array

            # Save the degrees of freedom.
            df2 = n_covariates * 2

            # Remove interaction column one by one.
            for cov_index in range(n_covariates):
                # Model the null matrix (without the interaction term).
                config_null, X_null, betas_null, rss_null = \
                    self.find_best_config(configs=null_model_configs,
                                          n_samples=n,
                                          n_covariates=n_covariates,
                                          cc_m=cc_m[:, mask],
                                          genotype=genotype,
                                          expression=y,
                                          exclude=cov_index
                                          )

                # Calculate the RSS just like Decon-eQTL does. This is not
                # how it was intended to work!
                # TODO delete this
                full_y = expr_m[row_index, :]
                null_y_hat = np.copy(alt_y_hat)
                null_y_hat[mask] = self.predict(X=X_null, betas=betas_null)
                rss_null = self.calc_rss(y=full_y, y_hat=null_y_hat)

                # Save the null model configuration.
                config_mask = [True if x != cov_index else False for x in range(n_covariates)]
                config_null_dict[cell_types_indices[cov_index]][row_index, config_mask] = config_null

                # Calculate the p-value.
                # TODO delete this
                p_value = self.calc_p_value(rss1=rss_null, rss2=rss_alt, df1=df2 - 1, df2=df2, n=full_y.shape[0])
                # p_value = self.calc_p_value(rss1=rss_null, rss2=rss_alt, df1=df2 - 1, df2=df2, n=n)

                if eqtl_indices[row_index] == "ENSG00000013573.17_12:31073901:rs7953706:T_A":
                    print(cell_types_indices[cov_index], rss_null, rss_alt, df2 - 1, df2, n, p_value)

                # Save the p-value.
                output_m[row_index, cov_index] = p_value

            # Print update for user.
            now_time = int(time.time())
            if last_print_time is None or (now_time - last_print_time) >= self.print_interval or (row_index + 1) == n_eqtls:
                print("\t[{}] {}/{} ieQTLs analysed [{:.2f}%]".format(time.strftime('%H:%M:%S', time.gmtime(now_time - start_time)),
                                                                      (row_index + 1),
                                                                      n_eqtls,
                                                                      (100 / n_eqtls) * (row_index + 1)))
                last_print_time = now_time
        print("")

        ########################################################################

        print("### STEP 4 ###")
        print("Saving results.")
        config_alt_df = pd.DataFrame(config_alt_m, index=eqtl_indices, columns=cell_types_indices)
        self.save_file(df=config_alt_df, outpath=os.path.join(self.outdir, "configuration_alt.txt.gz"))

        for cell_type, config_null_m in config_null_dict.items():
            config_null_df = pd.DataFrame(config_null_m, index=eqtl_indices, columns=cell_types_indices)
            self.save_file(df=config_null_df, outpath=os.path.join(self.outdir, "configuration_null_{}.txt.gz".format(cell_type)))

        output_df = pd.DataFrame(output_m,
                                 index=eqtl_indices,
                                 columns=["{}_pvalue".format(x) for x in cell_types_indices] +
                                         ["Beta{}_{}".format(i+1, x) for i, x in enumerate(cell_types_indices)] +
                                         ["Beta{}_{}:GT".format(len(cell_types_indices) + i + 1, x) for i, x in enumerate(cell_types_indices)])
        print(output_df)
        self.save_file(df=output_df, outpath=os.path.join(self.outdir, "deconvolutionResultsTest.txt.gz"))
        print("")

        ########################################################################

        pd.options.display.float_format = '{:.6f}'.format

        decon_df = self.load_file("../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv", header=0, index_col=0)
        for index in output_df.index:
            if index in decon_df.index:
                decon_result = decon_df.loc[[index], :].T
                decon_result.index = ["Beta1_CellMapNNLS_Astrocyte" if x == "Beta1_CellMapNNLS_Astrocyt" else x for x in decon_result.index]
                my_result = output_df.loc[[index], :].T
                comparison_df = decon_result.merge(my_result, left_index=True, right_index=True)
                comparison_df.columns = ["Decon-eQTL", "Martijn"]
                comparison_df["delta"] = (comparison_df["Decon-eQTL"] - comparison_df["Martijn"]).abs()
                print(comparison_df)
                exit()

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def validate_data(geno_df, expr_df, cc_df):
        if not geno_df.columns.tolist() == expr_df.columns.tolist():
            print("The genotype file header does not match the "
                  "expression file header.")
            exit()

        if not geno_df.columns.tolist() == cc_df.columns.tolist():
            print("The genotype file header does not match the "
                  "cell count file header.")
            exit()

    def find_best_config(self, configs, n_samples, n_covariates, cc_m,
                         genotype, expression, exclude=None):
        n_columns = n_covariates * 2
        if exclude is not None:
            n_columns -= 1

        # Create the model matrix. Leave the interaction columns blank.
        X = np.empty((n_samples, n_columns), np.float32)
        for cov_index in range(n_covariates):
            X[:, cov_index] = cc_m[cov_index, :]

        # Try different configurations for the genotype encoding.
        top_config = None
        top_X = None
        top_betas = None
        top_rss = np.inf
        for config in configs:
            # Fill in the alternative matrix with the right configuration
            # of allele encoding.
            for cov_index, flip in enumerate(config):
                cc_index = cov_index
                if exclude is not None and cc_index >= exclude:
                    cc_index += 1

                if flip:
                    X[:, n_covariates + cov_index] = (2 - genotype) * X[:, cc_index]
                else:
                    X[:, n_covariates + cov_index] = genotype * X[:, cc_index]

            # Calculate the R-squared.
            betas = self.fit(X, expression)
            # rss_alt = rnorm_alt * rnorm_alt
            rss_alt = self.calc_rss(y=expression, y_hat=self.predict(X=X, betas=betas))

            if rss_alt < top_rss:
                top_config = config
                top_X = np.copy(X)
                top_betas = betas
                top_rss = rss_alt

        return top_config, top_X, top_betas, top_rss

    @staticmethod
    def fit(X, y):
        betas, _ = nnls(X, y)
        return betas

    @staticmethod
    def predict(X, betas):
        return np.dot(X, betas)

    def fit_and_predict(self, X, y):
        return self.predict(X=X, betas=self.fit(X=X, y=y))

    @staticmethod
    def calc_pearsonr(x, y):
        x_dev = x - np.mean(x)
        y_dev = y - np.mean(y)
        dev_sum = np.sum(x_dev * y_dev)
        x_rss = np.sum(x_dev * x_dev)
        y_rss = np.sum(y_dev * y_dev)
        return dev_sum / np.sqrt(x_rss * y_rss)

    @staticmethod
    def calc_rss(y, y_hat):
        res = y - y_hat
        res_squared = res * res
        return np.sum(res_squared)

    @staticmethod
    def calc_p_value(rss1, rss2, df1, df2, n):
        """
        last row == stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2)) but this
        is faster somehow.

        https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html

        1 - I(a,b,x) = I(b, a, 1-x)
        """
        if rss2 >= rss1:
            return 1
        dfn = df2 - df1
        dfd = n - df2
        f_value = ((rss1 - rss2) / dfn) / (rss2 / dfd)
        p_value = betainc(dfd / 2, dfn / 2, 1 - ((dfn * f_value) / ((dfn * f_value) + dfd)))
        if p_value == 0:
            p_value = 2.2250738585072014e-308
        return p_value

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath), df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Missing genotype: {}".format(self.missing_geno))
        print("  > Alpha: {}".format(self.alpha))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
