#!/usr/bin/env python3

"""
File:         decon_eqtl_with_permutation_fdr.py
Created:      2021/07/12
Last Changed: 2021/07/14
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
import random
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
./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/expression_table/2020-07-16-MetaBrainDeconQtlGenes.TMM.SampSelect.ZeroVarRemov.covRemoved.expAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt.gz -stc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_cohort_matrix/sample_to_cohort.txt.gz -of CortexEUR-ciswoPermFDR-OLD -p 0

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt -stc /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/CortexEUR-cis/data/SampleToCohorts.txt.gz -of CortexEUR-cis-woPermFDR-New -p 0

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table_OLD.txt -stc /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/CortexEUR-cis/data/SampleToCohorts.txt.gz -of CortexEUR-cis-woPermFDR-New-OldCC -p 0
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
                  "of any kind. The methodology is based on Decon-eQTL: " \
                  "https://doi.org/10.1186/s12859-020-03576-5. Novel additions" \
                  "are the proper handling of missing values as well as " \
                  "calculating a permutation based FDR.".format(__program__,
                                                                __author__,
                                                                __license__)


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.stc_path = getattr(arguments, 'sample_to_cohort')
        self.nrows = getattr(arguments, 'rows')
        self.n_permutations = getattr(arguments, 'permutations')
        self.missing_geno = getattr(arguments, 'missing_genotype')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent.parent)
        self.outdir = os.path.join(outdir, "decon_eqtl_with_permutation_fdr", outfolder)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set other variables.
        self.alpha = 0.05
        self.print_interval = 30

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
                            help="The path to the genotype matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix.")
        parser.add_argument("-stc",
                            "--sample_to_cohort",
                            type=str,
                            required=True,
                            help="The path to the sample-to-cohort matri.x")
        parser.add_argument("-r",
                            "--rows",
                            type=int,
                            default=None,
                            help="The number of rows to analyze.")
        parser.add_argument("-p",
                            "--permutations",
                            type=int,
                            default=10,
                            help="The number of permutations to run.")
        parser.add_argument("-m",
                            "--missing_genotype",
                            type=int,
                            required=False,
                            default=-1,
                            help="The genotype value that equals a missing "
                                 "value. Default: -1. Note: has to be int.")
        parser.add_argument("-od",
                            "--outdir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the output path.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### STEP 1 ###")
        print("Loading data.")
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, nrows=self.nrows)
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, nrows=self.nrows)
        cc_df = self.load_file(self.cc_path, header=0, index_col=0)
        stc_df = self.load_file(self.stc_path, header=0, index_col=None)

        # Transpose if need be. We want samples always as columns.
        if cc_df.shape[0] == geno_df.shape[1] or cc_df.shape[0] == expr_df.shape[1]:
            print("\t  Transposing covariate matrix.")
            cc_df = cc_df.T

        # Remove columns with all nan. This line is just because
        # the expression file I used somehow had one column with nan's.
        expr_df.dropna(axis='columns', how='all', inplace=True)

        print("\tValidating input.")
        self.validate_data(geno_df=geno_df,
                           expr_df=expr_df,
                           cc_df=cc_df,
                           stc_df=stc_df)
        print("")

        ########################################################################

        print("### STEP 2 ###")
        print("Pre-processing data.")
        # Check if the requirements are met.
        if expr_df.min(axis=1).min() < 0:
            print("Error: expression matrix contains negative values.")
            exit()
        if cc_df.min(axis=1).min() < 0:
            print("Error: cell count matrix contains negative values.")
            exit()
        # TODO: check if the column sum is always == 1

        # Convert to numpy for speed.
        geno_m = geno_df.to_numpy(np.float64)
        expr_m = expr_df.to_numpy(np.float64)
        cc_m = cc_df.to_numpy(np.float64)
        
        # Save properties.
        eqtl_indices = expr_df.index + "_" + geno_df.index
        cell_types_indices = cc_df.index.to_numpy(dtype=object)
        del geno_df, expr_df, cc_df
        
        # Print info.
        n_samples = geno_m.shape[1]
        n_eqtls = geno_m.shape[0]
        n_covariates = cc_m.shape[0]
        print("\tN-samples: {}".format(n_samples))
        print("\tN-eQTLs: {}".format(n_eqtls))
        print("\tN-covariates: {}".format(n_covariates))
        print("")

        ########################################################################

        print("### STEP 3 ###")
        print("Creating permutation matrix.")
        perm_order_m = self.create_perm_orders(n_permutations=self.n_permutations,
                                               stc_df=stc_df)
        if self.n_permutations > 0:
            self.save_file(df=pd.DataFrame(perm_order_m, index=["permutation{}".format(i) for i in range(self.n_permutations)]),
                           outpath=os.path.join(self.outdir, "permutation_order.txt.gz"))
        n_permutation_values = n_eqtls * self.n_permutations
        print("\tN-permutation values in null distribution (per cell type): {}".format(n_permutation_values))
        print("")

        ########################################################################

        print("### STEP 4 ###")
        print("Analyzing.")

        # Initializing output matrices / arrays.
        real_pvalues_m = np.empty((n_eqtls, n_covariates), dtype=np.float64)
        perm_pvalues_m = np.empty((n_permutation_values, n_covariates), dtype=np.float64)
        betas_m = np.empty((n_eqtls, n_covariates * 2), dtype=np.float64)

        # Create a list of possible genotype encoding configuration. True means
        # we change the encoding by 2 - value. False means we do nothing.
        alt_model_configs = list(itertools.product([True, False], repeat=n_covariates))
        null_model_configs = list(itertools.product([True, False], repeat=n_covariates - 1))

        # Start loop.
        start_time = int(time.time())
        last_print_time = None
        for row_index in range(n_eqtls):
            # Print update for user.
            now_time = int(time.time())
            if last_print_time is None or (now_time - last_print_time) >= self.print_interval or row_index == (n_eqtls - 1):
                print("\t[{}] {}/{} ieQTLs analysed [{:.2f}%]".format(time.strftime('%H:%M:%S', time.gmtime(now_time - start_time)),
                                                                      row_index,
                                                                      (n_eqtls - 1),
                                                                      (100 / (n_eqtls - 1)) * row_index))
                last_print_time = now_time

            # Get the genotype.
            genotype = geno_m[row_index, :]

            # Construct the mask to remove missing values.
            mask = genotype != self.missing_geno
            n = np.sum(mask)

            # Model the alternative matrix (with the interaction term).
            # This is the matrix with expression ~ cc1 + cc2 + cc1 * geno +
            # cc2 * geno.
            config_alt, betas_alt, rss_alt = \
                self.model(
                    genotype=genotype[mask],
                    expression=expr_m[row_index, mask],
                    cell_fractions=cc_m[:, mask],
                    configs=alt_model_configs,
                    n_samples=n,
                    n_covariates=n_covariates
                )

            # Save the alternative model stats.
            # The beta's of the interaction terms are flipped if we
            # flipped the allele encoding. This makes it possible that some
            # betas are negative even though we use NNLS.
            flip_array = np.hstack((np.ones(n_covariates), np.vectorize({True: -1, False: 1}.get)(config_alt)))
            betas_m[row_index, :] = betas_alt * flip_array

            # Save the degrees of freedom of the alternative model.
            df2 = n_covariates * 2

            # Remove one interaction column (cc * geno) one by one and
            # determine the significance of the change in risdual sum of
            # squares with a f-test.
            for cov_index in range(n_covariates):

                # Model the null matrix (without the interaction term). In
                # this model 1 (and only 1!) of the cc * geno terms is removed.
                _, _, rss_null = \
                    self.model(
                        genotype=genotype[mask],
                        expression=expr_m[row_index, mask],
                        cell_fractions=cc_m[:, mask],
                        configs=null_model_configs,
                        n_samples=n,
                        n_covariates=n_covariates,
                        exclude=cov_index
                    )

                # Calculate and save the p-value.
                p_value = self.calc_p_value(rss1=rss_null, rss2=rss_alt, df1=df2 - 1, df2=df2, n=n)
                real_pvalues_m[row_index, cov_index] = p_value

                # Perform n permutations.
                for perm_index in range(self.n_permutations):
                    # Shuffle the genotype.
                    permutation_order = perm_order_m[perm_index, :]
                    shuffled_genotype = genotype[permutation_order][mask]

                    # Model the alternative matrix (with the interaction term)
                    # and shuffle the genotype of the interaction of interest.
                    _, _, rss_perm = \
                        self.model(
                            genotype=genotype[mask],
                            expression=expr_m[row_index, mask],
                            cell_fractions=cc_m[:, mask],
                            configs=alt_model_configs,
                            n_samples=n,
                            n_covariates=n_covariates,
                            shuffle_index=cov_index,
                            shuffled_genotype=shuffled_genotype,
                        )

                    # Calculate and save the permutation p-value.
                    perm_pvalue = self.calc_p_value(rss1=rss_null, rss2=rss_perm, df1=df2 - 1, df2=df2, n=n)
                    perm_pvalues_m[row_index * self.n_permutations + perm_index, cov_index] = perm_pvalue

        # Cap the p-values.
        real_pvalues_m[real_pvalues_m == 0] = 2.2250738585072014e-308
        print("")

        ########################################################################

        print("### STEP 5 ###")
        print("Calculating the permutation FDR.")
        perm_fdr_m = np.empty_like(real_pvalues_m)
        perm_fdr_m[:] = np.nan
        if self.n_permutations > 0:
            print("\nN-interaction (FDR < {}):".format(self.alpha))
            cov_length = np.max([len(x) for x in cell_types_indices])
            n_hits_total = 0
            for cov_index in range(n_covariates):
                n_hits = 0
                for row_index in range(n_eqtls):
                    # Get the real p-value.
                    real_pvalue = real_pvalues_m[row_index, cov_index]

                    # Get the rank of this p-value in both distributions.
                    rank = np.sum(real_pvalues_m[:, cov_index] <= real_pvalue)
                    perm_rank = np.sum(perm_pvalues_m[:, cov_index] <= real_pvalue)

                    # Calculate and safe the fdr.
                    fdr = (perm_rank / self.n_permutations) / rank
                    perm_fdr_m[row_index, cov_index] = fdr

                    # Increment the counter.
                    if fdr < self.alpha:
                        n_hits += 1

                print("\t{:{}s}  {}".format(cell_types_indices[cov_index], cov_length, n_hits))
                n_hits_total += n_hits

            print("\t{}".format("".join(["-"] * cov_length)))
            print("\t{:{}s}  {}".format("total", cov_length, n_hits_total))

            # Cap the permutation FDR values.
            perm_fdr_m[perm_fdr_m > 1] = 1
            perm_fdr_m[perm_fdr_m == 0] = 2.2250738585072014e-308
        print("")

        ########################################################################

        print(pd.DataFrame(real_pvalues_m, columns=["{}_pvalue".format(x) for x in cell_types_indices]))
        print(pd.DataFrame(perm_fdr_m, columns=["{}_FDR".format(x) for x in cell_types_indices]))
        print(pd.DataFrame(betas_m, columns=["Beta{}_{}".format(i+1, x) for i, x in enumerate(cell_types_indices)] + ["Beta{}_{}:GT".format(len(cell_types_indices) + i + 1, x) for i, x in enumerate(cell_types_indices)]))

        print("### STEP 6 ###")
        print("Saving results.")

        perm_pvalues_df = pd.DataFrame(perm_pvalues_m, columns=cell_types_indices)
        self.save_file(df=perm_pvalues_df, outpath=os.path.join(self.outdir, "permutation_pvalues.txt.gz"))

        output_df = pd.DataFrame(np.hstack((real_pvalues_m, perm_fdr_m, betas_m)),
                                 index=eqtl_indices,
                                 columns=["{}_pvalue".format(x) for x in cell_types_indices] +
                                         ["{}_FDR".format(x) for x in cell_types_indices] +
                                         ["Beta{}_{}".format(i+1, x) for i, x in enumerate(cell_types_indices)] +
                                         ["Beta{}_{}:GT".format(len(cell_types_indices) + i + 1, x) for i, x in enumerate(cell_types_indices)])
        print(output_df)
        self.save_file(df=output_df, outpath=os.path.join(self.outdir, "deconvolutionResults.txt.gz"))

        print("")

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
    def validate_data(geno_df, expr_df, cc_df, stc_df):
        if geno_df.columns.tolist() != expr_df.columns.tolist():
            print("The genotype file header does not match the "
                  "expression file header.")
            exit()

        if geno_df.columns.tolist() != cc_df.columns.tolist():
            print("The genotype file header does not match the "
                  "cell count file header.")
            exit()

        if geno_df.columns.tolist() != stc_df.iloc[:, 0].tolist():
            print("The sample-to-cohort file does not match the "
                  "genotype / expression file header.")
            exit()

    @staticmethod
    def create_perm_orders(n_permutations, stc_df):
        n_samples = stc_df.shape[0]
        cohorts = stc_df.iloc[:, 1].unique()
        valid_order = set([x for x in range(n_samples)])
        perm_order_m = np.empty((n_permutations, n_samples), dtype=np.uint8)
        for i in range(n_permutations):
            sample_indices = np.array([x for x in range(n_samples)])
            for cohort in cohorts:
                mask = (stc_df.iloc[:, 1] == cohort).to_numpy(dtype=bool)
                if np.sum(mask) <= 1:
                    continue
                copy = sample_indices[mask]
                random.shuffle(copy)
                sample_indices[mask] = copy

            if bool(valid_order.symmetric_difference(set(sample_indices))):
                print("Unvalid permutation order.")
                exit()

            perm_order_m[i, :] = sample_indices
        return perm_order_m

    def model(self, genotype, expression, cell_fractions, configs, n_samples,
              n_covariates, exclude=None, shuffle_index=None, shuffled_genotype=None):
        n_columns = n_covariates * 2
        if exclude is not None:
            n_columns -= 1

        # Create the model matrix. Leave the interaction columns blank for now.
        X = np.empty((n_samples, n_columns), np.float32)
        for cov_index in range(n_covariates):
            X[:, cov_index] = cell_fractions[cov_index, :]

        # Try different configurations for the genotype encoding.
        top_config = None
        top_betas = None
        top_rss = np.inf
        for config in configs:
            # Fill in the alternative matrix with the right configuration
            # of allele encoding.
            for cov_index, flip in enumerate(config):
                # If we exclude an interaction term we still have that cell
                # type fraction in the matrix. Therefore, when matching
                # interaction column position with cell fraction column
                # (cc_index) position we need to increment with 1 if the
                # cc_index with 1.
                cc_index = cov_index
                if exclude is not None and cc_index >= exclude:
                    cc_index += 1

                # Use a shuffled genotype vector if we are doing a permutation
                # analysis.
                genotype_a = genotype
                if cov_index == shuffle_index:
                    genotype_a = shuffled_genotype

                # Calculate genotype * cell fraction. If flip is true we
                # change the allele encoding (0 = 2, 1 = 1, 2 = 0).
                if flip:
                    X[:, n_covariates + cov_index] = (2 - genotype_a) * X[:, cc_index]
                else:
                    X[:, n_covariates + cov_index] = genotype_a * X[:, cc_index]

            # Model the expression vector as non-negative linear combination of
            # the model matrix. Determine the beta's as well as the squared
            # Euclidean norm. Note: squared Eucledian norm is identical to
            # calculating the RSS.
            betas, rnorm = self.fit(X, expression)
            rss_alt = rnorm * rnorm

            # Only safe the best configuration.
            if rss_alt < top_rss:
                top_config = config
                top_betas = betas
                top_rss = rss_alt

        return top_config, top_betas, top_rss

    @staticmethod
    def fit(X, y):
        return nnls(X, y)

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
        print("  > Sample-to-cohort path: {}".format(self.stc_path))
        print("  > N rows: {}".format(self.nrows))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > Missing genotype: {}".format(self.missing_geno))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
