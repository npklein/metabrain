#!/usr/bin/env python3

"""
File:         decon_eqtl_with_permutation_fdr.py
Created:      2021/07/12
Last Changed: 2021/08/05
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
import warnings
import argparse
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
./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/expression_table/2020-07-16-MetaBrainDeconQtlGenes.TMM.SampSelect.ZeroVarRemov.covRemoved.expAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_cohort_matrix/sample_to_dataset.txt.gz -of cortex_eur_cis

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/SampleToDataset.txt.gz -of CortexEUR-cis

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis-HalfNormalised/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceHalfNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-HalfNormalised/data/SampleToDataset.txt.gz -maf 5 -of CortexEUR-cis-HalfNormalizedMAF1

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/SampleToDataset.txt.gz -of CortexEUR-cis-PrimaryeQTLs -r 11803 -p 0

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR/genotype_table.txt -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/preprocess_scripts/pre_process_expression_matrix/CortexAFR-cis/data/SampleToDataset.txt.gz -of CortexAFR-cis-Replication-EUR

./decon_eqtl_with_permutation_fdr.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR-HalfNormalized/genotype_table.txt -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/decon_eqtl_replication_select_and_harmonize/CortexAFR-cis-Replication-EUR-HalfNormalized/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ForceHalfNormalised.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexAFR-cis-HalfNormalised/data/SampleToDataset.txt.gz -of CortexAFR-cis-Replication-EUR-HalfNormalised -maf 5 
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
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.nrows = getattr(arguments, 'rows')
        self.n_permutations = getattr(arguments, 'permutations')
        self.permutation_index_offset = getattr(arguments, 'permutation_index_offset')
        self.maf = getattr(arguments, 'minor_allele_frequency')
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
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-r",
                            "--rows",
                            type=int,
                            default=None,
                            help="The number of rows to analyze.")
        parser.add_argument("-p",
                            "--permutations",
                            type=int,
                            default=0,
                            help="The number of permutations to run. "
                                 "Default: 0.")
        parser.add_argument("-po",
                            "--permutation_index_offset",
                            type=int,
                            default=0,
                            help="The first index of the first permutation "
                                 "run. Default: 0.")
        parser.add_argument("-maf",
                            "--minor_allele_frequency",
                            type=int,
                            default=1,
                            help="The minimal required MAF for a genotype."
                                 "Default: 1%")
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
        std_df = self.load_file(self.std_path, header=0, index_col=None)

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
                           std_df=std_df)
        print("", flush=True)

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
        std_m = std_df.to_numpy(object)

        # Replace missing values with nan
        geno_m[geno_m == self.missing_geno] = np.nan

        # Save properties.
        eqtl_indices = expr_df.index + "_" + geno_df.index
        cell_types_indices = cc_df.index.to_numpy(dtype=object)
        dataset_counts = list(zip(*np.unique(std_m[:, 1], return_counts=True)))
        dataset_counts.sort(key=lambda x: -x[1])
        datasets = [x[0] for x in dataset_counts]
        dataset_sample_sizes = np.array([x[1] for x in dataset_counts])
        if np.min(dataset_sample_sizes) <= 1:
            print("Dataset with 0 or 1 sample is not allowed.")
            exit()
        del geno_df, expr_df, cc_df, std_df

        # Calculate MAF and filter on cut-off.
        maf_a = self.calculate_maf(m=geno_m)
        maf_mask = maf_a > (self.maf / 100)
        if np.sum(maf_mask) == 0:
            print("No genotype above MAF cut-off")
            exit()
        if np.sum(maf_mask) != np.size(maf_mask):
            eqtl_indices = eqtl_indices[maf_mask]
            geno_m = geno_m[maf_mask, :]
            expr_m = expr_m[maf_mask, :]

        # Calculate missing values stats. Set 1 as np.nan for a better average.
        prcnt_missing_m = np.sum(np.isnan(geno_m), axis=1) / geno_m.shape[1]
        prcnt_missing_m[prcnt_missing_m == 1] = np.nan

        prcnt_missing_per_dataset_m = np.empty((geno_m.shape[0], len(datasets)), dtype=np.float64)
        for dataset_index, dataset in enumerate(datasets):
            dataset_mask = std_m[:, 1] == dataset
            prcnt_missing_per_dataset_m[:, dataset_index] = np.sum(np.isnan(geno_m[:, dataset_mask]), axis=1) / np.sum(dataset_mask)
        prcnt_missing_per_dataset_m[prcnt_missing_per_dataset_m == 1] = np.nan

        # Print info.
        n_eqtls = geno_m.shape[0]
        n_samples = geno_m.shape[1]
        n_covariates = cc_m.shape[0]
        print("\tN-eQTLs: {:,}".format(n_eqtls))
        print("\tN-samples: {:,}".format(n_samples))
        print("\tN-covariates: {:,}".format(n_covariates))
        print("\tN-datasets: {:,}".format(len(datasets)))
        if np.sum(maf_mask) != np.size(maf_mask):
            print("\tMAF% (before filter): N = {:,}\tmean = {:.2f}\tmin = {:.2f}\tmax = {:.2f}".format(np.size(maf_a), np.mean(maf_a), np.min(maf_a), np.max(maf_a)))
            print("\tMAF% (after filter):  N = {:,}\tmean = {:.2f}\tmin = {:.2f}\tmax = {:.2f}".format(np.size(maf_a[maf_mask]), np.mean(maf_a[maf_mask]), np.min(maf_a[maf_mask]), np.max(maf_a[maf_mask])))
        else:
            print("\tMAF%: mean = {:.2f}\tmin = {:.2f}\tmax = {:.2f}".format(np.mean(maf_a), np.min(maf_a), np.max(maf_a)))
        print("\t%-missing (avg per eQTL): {:.2f}%".format(np.nanmean(prcnt_missing_m) * 100))
        print("\t%-missing (per eQTL per dataset): mean = {:.2f}%\tmin = {:.2f}\tmax = {:.2f}%".format(np.mean(np.nanmean(prcnt_missing_per_dataset_m, axis=1)) * 100, np.nanmin(prcnt_missing_per_dataset_m) * 100, np.nanmax(prcnt_missing_per_dataset_m) * 100))
        print("", flush=True)

        ########################################################################

        print("### STEP 3 ###")
        print("Analyzing eQTLs.")

        # Calculate and print some info about the analyses to be performed.
        n_configurations_alt = (n_covariates * 2) + 2
        n_configurations_null = ((n_covariates - 1) * 2) + 2
        n_models = n_eqtls * (n_configurations_alt + n_covariates * n_configurations_null)
        print("\tN-configurations (full model): {:,}".format(n_configurations_alt))
        print("\tN-configurations (cell type model): {:,}".format(n_configurations_null))
        print("\tN-models: {:,}".format(n_models))
        print("")

        # Initializing output matrices / arrays.
        real_pvalues_m = np.empty((n_eqtls, n_covariates), dtype=np.float64)
        betas_m = np.empty((n_eqtls, n_covariates * 2), dtype=np.float64)
        rss_null_m = np.empty((n_eqtls, n_covariates), dtype=np.float64)

        # Create a list of possible genotype encoding configuration. True means
        # we change the encoding by 2 - value. False means we do nothing.
        alt_model_configs = self.create_model_configs(n=n_covariates)
        null_model_configs = self.create_model_configs(n=n_covariates - 1)

        # Save the degrees of freedom the alternative model.
        df = n_covariates * 2

        # Start loop.
        start_time = int(time.time())
        last_print_time = None
        for row_index in range(n_eqtls):
            # Print update for user.
            now_time = int(time.time())
            if n_eqtls > 1 and (last_print_time is None or (now_time - last_print_time) >= self.print_interval or row_index == (n_eqtls - 1)):
                print("\t[{}] {}/{} eQTLs analysed [{:.2f}%]".format(time.strftime('%H:%M:%S', time.gmtime(now_time - start_time)),
                                                                     row_index,
                                                                     (n_eqtls - 1),
                                                                     (100 / (n_eqtls - 1)) * row_index),
                      flush=True)
                last_print_time = now_time

            # Get the genotype.
            genotype = geno_m[row_index, :]

            # Construct the mask to remove missing values.
            mask = ~np.isnan(genotype)
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

            # Remove one interaction column (cc * geno) one by one and
            # determine the significance of the change in residuals sum of
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
                p_value = self.calc_p_value(rss1=rss_null,
                                            rss2=rss_alt,
                                            df1=df - 1,
                                            df2=df,
                                            n=n)
                real_pvalues_m[row_index, cov_index] = p_value

                # Save the RSS null for permutations later.
                rss_null_m[row_index, cov_index] = rss_null

        # Cap the p-values.
        real_pvalues_m[real_pvalues_m == 0] = 2.2250738585072014e-308

        # Print the number of significant hits.
        self.print_n_signif(m=real_pvalues_m, colnames=cell_types_indices, type="p-value")
        print("", flush=True)

        ########################################################################

        if self.permutation_index_offset == 0:
            print("### STEP 4 ###")
            print("Saving results.")

            self.save_file(df=pd.DataFrame(maf_a[maf_mask], index=eqtl_indices, columns=["NAF"]), outpath=os.path.join(self.outdir, "MAF.txt.gz"))

            # Combine missing matrices and set np.nan == 1 back.
            combined_prcnt_missing_m = np.hstack((prcnt_missing_m[:, np.newaxis], prcnt_missing_per_dataset_m))
            combined_prcnt_missing_m[np.isnan(combined_prcnt_missing_m)] = 1
            combined_prcnt_missing_df = pd.DataFrame(combined_prcnt_missing_m, columns=["Total"] + datasets, index=eqtl_indices)
            print(combined_prcnt_missing_df)
            self.save_file(df=combined_prcnt_missing_df, outpath=os.path.join(self.outdir, "missing_values_info.txt.gz"))

            decon_df = pd.DataFrame(np.hstack((real_pvalues_m, betas_m)),
                                    index=eqtl_indices,
                                    columns=["{}_pvalue".format(x) for x in cell_types_indices] +
                                            ["Beta{}_{}".format(i+1, x) for i, x in enumerate(cell_types_indices)] +
                                            ["Beta{}_{}:GT".format(len(cell_types_indices) + i + 1, x) for i, x in enumerate(cell_types_indices)])
            print(decon_df)
            self.save_file(df=decon_df, outpath=os.path.join(self.outdir, "deconvolutionResults.txt.gz"))
            del prcnt_missing_m, prcnt_missing_per_dataset_m, combined_prcnt_missing_df, decon_df, betas_m

            print("", flush=True)

        #######################################################################

        if self.n_permutations <= 0:
            exit()

        print("### STEP 5 ###")
        print("Performing permutations.")

        # Calculate and print some info about the analyses to be performed.
        n_permutation_values = n_eqtls * self.n_permutations
        n_perm_models = n_eqtls * (self.n_permutations * n_covariates * n_configurations_null)
        print("\tN-permutation p-values (per cell type): {:,}".format(n_permutation_values))
        print("\tN-models: {:,}".format(n_perm_models))
        print("")

        # Create the permutation orders.
        perm_order_m = self.create_perm_orders(n_permutations=self.n_permutations,
                                               n_samples=n_samples,
                                               datasets=datasets,
                                               std_m=std_m,
                                               seed_offset=self.permutation_index_offset)

        # Calculate the overlap between the genotype matrices
        perm_overlap_m = self.check_permutation_overlap(geno_m=geno_m, perm_order_m=perm_order_m)

        geno_nan_mask = np.isnan(geno_m)
        nanfilled_geno_m = np.copy(geno_m)
        if np.sum(geno_nan_mask) > 0:
            # Calculate the dataset genotype means per eQTL and fill nan values
            # with the dataset mean.
            geno_dataset_mean_m = self.calculate_geno_mean_per_dataset(geno_m=geno_m,
                                                                       datasets=datasets,
                                                                       std_m=std_m)
            nanfilled_geno_m[geno_nan_mask] = geno_dataset_mean_m[geno_nan_mask]

        # Initializing output matrices / arrays.
        perm_pvalues_m = np.empty((n_eqtls, n_covariates, self.n_permutations), dtype=np.float64)

        # Start loop.
        start_time = int(time.time())
        last_print_time = None
        for row_index in range(n_eqtls):
            # Print update for user.
            now_time = int(time.time())
            if n_eqtls > 1 and (last_print_time is None or (now_time - last_print_time) >= self.print_interval or row_index == (n_eqtls - 1)):
                print("\t[{}] {}/{} eQTLs analysed [{:.2f}%]".format(time.strftime('%H:%M:%S', time.gmtime(now_time - start_time)),
                                                                     row_index,
                                                                     (n_eqtls - 1),
                                                                     (100 / (n_eqtls - 1)) * row_index),
                      flush=True)
                last_print_time = now_time

            # Get the genotype arrays.
            genotype = geno_m[row_index, :]
            nanfilled_genotype = nanfilled_geno_m[row_index, :]

            # Construct the mask to remove missing values.
            mask = ~np.isnan(genotype)
            n = np.sum(mask)

            # Loop over the covariates.
            for cov_index in range(n_covariates):

                # Perform n permutations.
                for perm_index in range(self.n_permutations):
                    # Calculate the interaction term for both allele encodings.
                    # Use the NaN-filled genotype array for this.
                    interaction_a = nanfilled_genotype * cc_m[cov_index, :]
                    interaction_flipped_a = (2 - nanfilled_genotype) * cc_m[cov_index, :]

                    # Shuffle the indices.
                    perm_order = perm_order_m[perm_index, :]
                    interaction_a = interaction_a[perm_order]
                    interaction_flipped_a = interaction_flipped_a[perm_order]

                    # Model the alternative matrix (with the interaction
                    # term) and shuffle the genotype of the interaction of
                    # interest.
                    _, _, rss_perm = \
                        self.model(
                            genotype=genotype[mask],
                            expression=expr_m[row_index, mask],
                            cell_fractions=cc_m[:, mask],
                            configs=alt_model_configs,
                            n_samples=n,
                            n_covariates=n_covariates,
                            shuffle_index=cov_index,
                            shuffle_inter=interaction_a[mask],
                            shuffle_inter_flipped=interaction_flipped_a[mask]
                        )

                    # Calculate and save the permutation p-value.
                    perm_pvalue = self.calc_p_value(
                        rss1=rss_null_m[row_index, cov_index],
                        rss2=rss_perm,
                        df1=df - 1,
                        df2=df,
                        n=n)
                    perm_pvalues_m[row_index, cov_index, perm_index] = perm_pvalue

        print("", flush=True)

        # #######################################################################
        #
        # print("### STEP 6 ###")
        # print("Calculating permutation based FDR.")
        # perm_fdr_m = np.empty_like(real_pvalues_m)
        # perm_fdr_m[:] = np.nan
        # for cov_index in range(n_covariates):
        #     for row_index in range(n_eqtls):
        #         # Get the real p-value.
        #         real_pvalue = real_pvalues_m[row_index, cov_index]
        #
        #         # Get the rank of this p-value in both distributions.
        #         rank = np.sum(real_pvalues_m[:, cov_index] <= real_pvalue)
        #         perm_rank = np.sum(perm_pvalues_m[:, cov_index, :] <= real_pvalue)
        #
        #         # Calculate and safe the fdr.
        #         perm_fdr_m[row_index, cov_index] = (perm_rank / self.n_permutations) / rank
        #
        # # Cap the permutation FDR values.
        # perm_fdr_m[perm_fdr_m > 1] = 1
        # perm_fdr_m[perm_fdr_m == 0] = 2.2250738585072014e-308
        #
        # # Print the number of significant hits.
        # self.print_n_signif(m=perm_fdr_m, colnames=cell_types_indices, type="fdr-value")
        # print("", flush=True)
        #
        # ########################################################################

        print("### STEP 7 ###")
        print("Saving results.")

        self.save_matrix(m=perm_pvalues_m, outpath=os.path.join(self.outdir, "permutation_pvalues_{}_until_{}.npy".format(self.permutation_index_offset, self.permutation_index_offset + self.n_permutations - 1)))
        self.save_matrix(m=perm_order_m, outpath=os.path.join(self.outdir, "perm_orders_{}_until_{}.npy".format(self.permutation_index_offset, self.permutation_index_offset + self.n_permutations - 1)))
        self.save_matrix(m=perm_overlap_m, outpath=os.path.join(self.outdir, "perm_order_overlap_{}_until_{}.npy".format(self.permutation_index_offset, self.permutation_index_offset + self.n_permutations - 1)))

        lowest_pvalues_m = np.transpose(np.min(perm_pvalues_m, axis=0))
        lowest_pvalues_df = pd.DataFrame(lowest_pvalues_m,
                                         index=["permutation_{:06d}".format(self.permutation_index_offset + perm_index) for perm_index in range(self.n_permutations)],
                                         columns=["{}_pvalue".format(cell_type) for cell_type in cell_types_indices])
        self.save_file(df=lowest_pvalues_df, outpath=os.path.join(self.outdir, "lowest_permutation_pvalues_{}_until_{}.txt.gz".format(self.permutation_index_offset, self.permutation_index_offset + self.n_permutations - 1)))

        # perm_fdr_df = pd.DataFrame(perm_fdr_m, columns=["{}_FDR".format(cell_type) for cell_type in cell_types_indices])
        # self.save_file(df=perm_fdr_df, outpath=os.path.join(self.outdir, "permutation_FDR.txt.gz"))

        print("", flush=True)

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
    def validate_data(geno_df, expr_df, cc_df, std_df):
        if geno_df.columns.tolist() != expr_df.columns.tolist():
            print("The genotype file header does not match the "
                  "expression file header.")
            exit()

        if geno_df.columns.tolist() != cc_df.columns.tolist():
            print("The genotype file header does not match the "
                  "cell count file header.")
            exit()

        if geno_df.columns.tolist() != std_df.iloc[:, 0].tolist():
            print("The sample-to-dataset file does not match the "
                  "genotype / expression file header.")
            exit()

    @staticmethod
    def calculate_maf(m):
        rounded_m = np.copy(m)
        rounded_m = np.rint(rounded_m)

        zero_a = np.sum(rounded_m == 0, axis=1)
        one_a = np.sum(rounded_m == 1, axis=1)
        two_a = np.sum(rounded_m == 2, axis=1)

        allele1_a = (zero_a * 2) + one_a
        allele2_a = (two_a * 2) + one_a

        return np.minimum(allele1_a, allele2_a) / (allele1_a + allele2_a)

    @staticmethod
    def create_model_configs(n):
        """
        Create the allele encoding configurations. All configurations could
        be created using list(itertools.product([True, False], repeat=n))
        however they mention in the article that Decon-eQTL restrict the
        configurations to max one opposite. This limits the
        configurations from k^2 to (2*k) + 2.

        Example for n = 3
        Output:
            [ False False False ]
            [ True False False ]
            [ False True False ]
            [ False False True ]
            [ True True True ]
            [ False True True ]
            [ True False True ]
            [ True True False ]

        """
        configurations = []

        false_array = np.zeros(n, dtype=bool)
        true_array = np.ones(n, dtype=bool)
        for value, array in zip([True, False], [false_array, true_array]):
            configurations.append(array)
            for i in range(n):
                configuration = np.copy(array)
                configuration[i] = value
                configurations.append(configuration)

        return configurations

    @staticmethod
    def create_perm_orders(n_permutations, n_samples, datasets, std_m, seed_offset):
        """
        Shuffles an array of size n_samples into a random order (with seed)
        N times. However, this function only shuffles samples within the same
        dataset.
        """
        default_order = np.array([x for x in range(n_samples)], dtype=np.uint16)
        all_indices = set(default_order)
        perm_order_m = np.empty((n_permutations, n_samples), dtype=np.uint16)
        for i in range(perm_order_m.shape[0]):
            perm_order = np.copy(default_order)
            for dataset in datasets:
                dataset_mask = std_m[:, 1] == dataset
                if np.sum(dataset_mask) <= 1:
                    continue

                copy = perm_order[dataset_mask]
                random.Random(seed_offset + i).shuffle(copy)
                perm_order[dataset_mask] = copy

            if set(perm_order) != all_indices:
                print("Invalid permutation order.")
                exit()

            perm_order_m[i, :] = perm_order
        return perm_order_m

    @staticmethod
    def check_permutation_overlap(geno_m, perm_order_m):
        """
        Function to calculate how identical the permuted genotype matrix is to
        the original genotype matrix. Returns a matrix with n-eQTLs x n-permutations
        with the % overlap in the cells.
        """
        perm_overlap_m = np.empty((geno_m.shape[0], perm_order_m.shape[0]), dtype=np.float64)

        for perm_index in range(perm_order_m.shape[0]):
            perm_order = perm_order_m[perm_index, :]
            perm_geno_m = np.copy(geno_m)
            perm_geno_m = perm_geno_m[:, perm_order]
            perm_overlap_m[:, perm_index] = np.sum(np.logical_and(geno_m == perm_geno_m, ~np.isnan(geno_m)), axis=1) / geno_m.shape[1]

        return perm_overlap_m

    @staticmethod
    def calculate_geno_mean_per_dataset(geno_m, datasets, std_m):
        """
        Method for calculating the mean genotype per dataset per eQTL. Missing
        values are not included. Returns a matrix with n-eQTLs x n-samples
        with the mean genotype per dataset in the cells.
        """
        geno_dataset_mean_m = np.empty_like(geno_m)

        # Ignore RuntimeWarning: Mean of empty slice
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            # Fill in the missing values with the dataset mean.
            for dataset_index, dataset in enumerate(datasets):
                dataset_mask = std_m[:, 1] == dataset
                geno_dataset_m = geno_m[:, dataset_mask]
                dataset_mean_a = np.nanmean(geno_dataset_m, axis=1)
                geno_dataset_mean_m[:, dataset_mask] = np.tile(dataset_mean_a[:, np.newaxis], np.sum(dataset_mask))

        return geno_dataset_mean_m

    @staticmethod
    def model(genotype, expression, cell_fractions, configs, n_samples,
              n_covariates, exclude=None, shuffle_index=None,
              shuffle_inter=None, shuffle_inter_flipped=None):
        """
        Create the interaction model. Try different allele encodings and
        find the optimal configurations. Only the best configuration is stored
        and returned.
        """
        if shuffle_index is not None and (shuffle_inter is None or shuffle_inter_flipped is None):
            print("Both shuffle_inter and shuffle_inter_flipped arguments "
                  "are required if shuffle index is set.")
            exit()

        n_columns = n_covariates * 2
        if exclude is not None:
            n_columns -= 1

        # Create the model matrix. Leave the interaction columns blank for now.
        X = np.empty((n_samples, n_columns), dtype=np.float64)
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
                # (cc_index) position we need to increment with 1 for all
                # cell fractions > exclude column.
                cc_index = cov_index
                if exclude is not None and cc_index >= exclude:
                    cc_index += 1

                # Use a shuffled genotype vector if we are doing a permutation
                # analysis.
                if cov_index == shuffle_index:
                    if flip:
                        X[:, n_covariates + cov_index] = shuffle_inter_flipped
                    else:
                        X[:, n_covariates + cov_index] = shuffle_inter
                else:
                    # Calculate genotype * cell fraction. If flip is true we
                    # change the allele encoding (0 = 2, 1 = 1, 2 = 0).
                    if flip:
                        X[:, n_covariates + cov_index] = (2 - genotype) * X[:, cc_index]
                    else:
                        X[:, n_covariates + cov_index] = genotype * X[:, cc_index]

            # Check if all values are positive.
            if np.min(X) < 0:
                print("Error: negative values in regression matrix.")
                exit()

            # Model the expression vector as non-negative linear combination of
            # the model matrix. Determine the beta's as well as the squared
            # Euclidean norm. Note: squared Eucledian norm is identical to
            # calculating the RSS.
            betas, rnorm = nnls(X, expression)
            rss_alt = rnorm * rnorm

            # Only safe the best configuration.
            if rss_alt < top_rss:
                top_config = config
                top_betas = betas
                top_rss = rss_alt

        return top_config, top_betas, top_rss

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
    def calc_pearsonr(x, y):
        x_dev = x - np.mean(x)
        y_dev = y - np.mean(y)
        dev_sum = np.sum(x_dev * y_dev)
        x_rss = np.sum(x_dev * x_dev)
        y_rss = np.sum(y_dev * y_dev)
        return dev_sum / np.sqrt(x_rss * y_rss)

    def print_n_signif(self, m, colnames, type):
        print("\nN-interaction ({} < {}):".format(type, self.alpha))
        n_hits_a = (m < self.alpha).sum(axis=0)
        n_hits_total = np.sum(n_hits_a)
        cov_length = np.max([len(x) for x in colnames])
        hits_length = np.max([len(str(x)) for x in n_hits_a] + [len(str(n_hits_total))])
        for n_hits, cell_type in zip(n_hits_a, colnames):
            print("\t{:{}s}  {:{}d}".format(cell_type, cov_length, n_hits, hits_length))
        print("\t{}".format("".join(["-"] * cov_length)))
        print("\t{:{}s}  {:{}d}".format("total", cov_length, n_hits_total, hits_length))

        print("", flush=True)

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath), df.shape))

    @staticmethod
    def save_matrix(m, outpath):
        with open(outpath, 'wb') as f:
            np.save(f, m)
        f.close()

        print("\tSaved matrix: {} "
              "with shape: {}".format(os.path.basename(outpath), m.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N rows: {}".format(self.nrows))
        print("  > N permutations: {}".format(self.n_permutations))
        print("  > Permutation index offset: {}".format(self.permutation_index_offset))
        print("  > MAF: {}%".format(self.maf))
        print("  > Missing genotype: {}".format(self.missing_geno))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
