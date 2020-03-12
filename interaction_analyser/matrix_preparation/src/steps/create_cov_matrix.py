"""
File:         create_cov_matrices.py
Created:      2020/03/12
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

# Local application imports.
from src.utilities import prepare_output_dir, check_file_exists
from src.df_utilities import load_dataframe, save_dataframe


class CreateCovMatrix:
    def __init__(self, settings, marker_file, force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param marker_file: string, path to the marker file.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.cov_file = settings["covariate_datafile"]
        self.tech_covs = settings["technical_covariates"]
        self.cohorts = settings["cohorts"]
        self.pheno_file = settings["phenotype_datafile"]
        self.eig_file = settings["eigenvectors_datafile"]
        self.n_eigen = settings["num_eigenvectors"]
        self.eig_bef_cov_corr_file = settings["eigenvectors_before_cov_corr_datafile"]
        self.marker_file = marker_file
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_cov_matrix')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "covariates-cortex.txt.gz")

        # Variables.
        self.covariates = None
        self.sex_dict = {"M": 0, "F": 1, np.nan: -1}

    def start(self):
        print("Starting creating covariate file.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.covariates = load_dataframe(inpath=self.outpath, header=0,
                                             index_col=0)
        else:
            self.covariates = self.combine_files()
            self.save()

    def combine_files(self):
        # read the covariates file.
        print("Loading covariate matrix.")
        cov_df = load_dataframe(inpath=self.cov_file, header=0, index_col=0)

        # filter technical covariates and the cohorts out.
        tech_cov_df = cov_df.loc[:, self.tech_covs]
        cohorts_df = cov_df.loc[:, self.cohorts]
        del cov_df

        # validate the cohorts.
        print("Validating cohorts.")
        colsums = cohorts_df.sum(axis=1)
        cohorts_df['ENA-EU'] = 0
        cohorts_df.loc[colsums == 0, 'ENA-EU'] = 1
        if not cohorts_df.sum(axis=1).all():
            print("\tSome samples do not have a cohort.")
            exit()
        else:
            print("\tValid.")

        # read the phenotype file.
        print("Loading phenotype matrix.")
        pheno_df = load_dataframe(inpath=self.pheno_file, header=0,
                                  index_col=4, low_memory=False)
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

        # read the marker genes expression file.
        print("Loading marker genes matrix.")
        marker_df = load_dataframe(self.marker_file, header=0, index_col=0)
        marker_df = marker_df.T

        # merge.
        print("Merging matrices.")
        tmp_df1 = tech_cov_df.merge(cohorts_df, left_index=True,
                                    right_index=True)
        tmp_df2 = tmp_df1.merge(gender_df, left_index=True, right_index=True)
        tmp_df3 = tmp_df2.merge(eigen_df, left_index=True, right_index=True)
        tmp_df4 = tmp_df3.merge(cov_cor_df, left_index=True, right_index=True)
        merged_df = tmp_df4.merge(marker_df, left_index=True, right_index=True)
        del tmp_df1, tmp_df2, tmp_df3, tmp_df4
        print("\tShape: {}".format(merged_df.shape))

        # post process the matrix.
        print("Post-process covariate matrix.")
        merged_df = merged_df.T
        merged_df.index = merged_df.index.set_names(['Sample'])
        merged_df.index.name = "-"
        print("\tShape: {}".format(merged_df.shape))

        return merged_df

    def save(self):
        save_dataframe(df=self.covariates, outpath=self.outpath,
                       index=True, header=True)

    def get_covariates(self):
        return self.covariates

    def get_sex_dict(self):
        return self.sex_dict

    def get_outpath(self):
        return self.outpath

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariates input file: {}".format(self.cov_file))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > Cohorts: {}".format(self.cohorts))
        print("  > Phenotype input file: {}".format(self.pheno_file))
        print("  > Eigenvectors input file: {}".format(self.eig_file))
        print("  > N. Eigenvectors: {}".format(self.n_eigen))
        print("  > Eigenvec before cov. corr. input file: {}".format(self.eig_bef_cov_corr_file))
        print("  > Markers input file: {}".format(self.marker_file))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
