"""
File:         create_cov_matrices.py
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
from functools import reduce
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.
from matrix_preparation.src.utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class CreateCovMatrix:
    def __init__(self, settings, cohort_df, decon_df, sample_order, force,
                 outdir):
        self.cov_file = settings["covariate_datafile"]
        self.tech_covs = settings["technical_covariates"]
        self.pheno_file = settings["phenotype_datafile"]
        self.eig_file = settings["eigenvectors_datafile"]
        self.n_eigen = settings["num_eigenvectors"]
        self.eig_bef_cov_corr_file = settings["eigenvectors_before_cov_corr_datafile"]
        self.cohort_df = cohort_df
        self.decon_df = decon_df
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_cov_matrix')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "covariates_table.txt.gz")

        # Variables.
        self.cov_df = None
        self.sex_dict = {"M": 0, "F": 1, np.nan: -1}

    def start(self):
        print("Starting creating covariate file.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or  self.force:
            self.cov_df = self.combine_files()
            self.save()

    def combine_files(self):
        # read the covariates file.
        print("Loading covariate matrix.")
        cov_df = load_dataframe(inpath=self.cov_file, header=0, index_col=0)
        tech_cov_df = cov_df[self.tech_covs].copy()
        del cov_df

        # read the phenotype file.
        print("Loading phenotype matrix.")
        pheno_df = load_dataframe(inpath=self.pheno_file, header=0,
                                  index_col=4, low_memory=False)

        # Combine the two gender columns, keep 'sex.by.expression' as main
        # gender ans use 'Gender' when no information is available.
        pheno_df = pheno_df.loc[:, ["Gender", "sex.by.expression"]]
        pheno_df.replace("no expression available", np.nan, inplace=True)
        pheno_df["SEX"] = pheno_df['sex.by.expression'].combine_first(
            pheno_df['Gender'])
        gender_df = pheno_df["SEX"].to_frame()
        del pheno_df
        gender_df = gender_df.replace({"SEX": self.sex_dict})

        # read the eigenvectors file.
        print("Loading eigenvectors matrix.")
        eigen_df = load_dataframe(self.eig_file, header=0, index_col=0)
        eigen_df = eigen_df.loc[:, ["Comp{}".format(x) for x in range(1, self.n_eigen + 1)]]

        # read the eigenvectors before covariate correction file.
        print("Loading eigenvectors before cov. correction matrix.")
        cov_cor_df = load_dataframe(self.eig_bef_cov_corr_file, header=0, index_col=0)
        cov_cor_df.columns = ["PC1-before-cov-correction", "PC2-before-cov-correction"]

        # merge.
        print("Merging matrices.")
        comb_cov = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
                          [tech_cov_df, self.cohort_df, gender_df, eigen_df, cov_cor_df, self.decon_df])
        comb_cov = comb_cov.T
        comb_cov = comb_cov[self.sample_order]
        comb_cov.index.name = "-"
        print("\tShape: {}".format(comb_cov.shape))

        # Remove old dataframes.
        del tech_cov_df, gender_df, eigen_df, cov_cor_df

        return comb_cov

    def save(self):
        save_dataframe(df=self.cov_df, outpath=self.outpath,
                       index=True, header=True)

    def clear_variables(self):
        self.cohorts_df = None
        self.decon_df = None
        self.sample_order = None
        self.force = None

    def print_arguments(self):
        print("Arguments:")
        print("  > Cohorts: {}".format(self.cohorts_df.shape))
        print("  > Deconvolution: {}".format(self.decon_df.shape))
        print("  > Covariates input file: {}".format(self.cov_file))
        print("  > Technical Covarates: {}".format(self.tech_covs))
        print("  > Phenotype input file: {}".format(self.pheno_file))
        print("  > Eigenvectors input file: {}".format(self.eig_file))
        print("  > N. Eigenvectors: {}".format(self.n_eigen))
        print("  > Eigenvec before cov. corr. input file: {}".format(self.eig_bef_cov_corr_file))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
