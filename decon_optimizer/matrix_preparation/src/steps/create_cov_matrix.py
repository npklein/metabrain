"""
File:         create_cov_matrices.py
Created:      2020/10/08
Last Changed: 2020/10/13
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
    def __init__(self, settings, tech_covs_file, tech_covs_df, cohort_file,
                 cohort_df, decon_file, decon_df, sample_dict, sample_order,
                 force, outdir):
        self.eig_file = settings["eigenvectors_datafile"]
        self.n_eigen = settings["num_eigenvectors"]
        self.tech_covs_file = tech_covs_file
        self.tech_covs_df = tech_covs_df
        self.cohort_file = cohort_file
        self.cohort_df = cohort_df
        self.decon_file = decon_file
        self.decon_df = decon_df
        self.sample_dict = sample_dict
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
        if not check_file_exists(self.outpath) or self.force:
            self.cov_df = self.combine_files()
            self.save()

    def combine_files(self):
        # loading technical covariates matrix.
        print("Loading technical covariates matrix.")
        if self.tech_covs_df is None:
            self.tech_covs_df = load_dataframe(self.tech_covs_file,
                                               header=0,
                                               index_col=0)

        # loading cohort matrix.
        print("Loading cohort matrix.")
        if self.cohort_df is None:
            self.cohort_df = load_dataframe(self.cohort_file,
                                            header=0,
                                            index_col=0)

        # read the eigenvectors file.
        print("Loading eigenvectors matrix.")
        eigen_df = load_dataframe(self.eig_file, header=0, index_col=0)
        eigen_df = eigen_df.loc[:, ["Comp{}".format(x) for x in range(1, self.n_eigen + 1)]]
        eigen_df.index = [self.sample_dict[x] if x in self.sample_dict else x for x in eigen_df.index]
        eigen_df = eigen_df.loc[self.sample_order, :]

        # loading deconvolution matrix.
        print("Loading deconvolution matrix.")
        if self.decon_df is None:
            self.decon_df = load_dataframe(self.decon_file,
                                           header=0,
                                           index_col=0)

        # merge.
        print("Merging matrices.")
        comb_cov = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
                          [self.tech_covs_df, self.cohort_df.T, eigen_df, self.decon_df])
        comb_cov = comb_cov.T
        comb_cov.index.name = "-"
        print("\tShape: {}".format(comb_cov.shape))

        # Validate sample order.
        if not comb_cov.columns.equals(self.sample_order):
            comb_cov = comb_cov[self.sample_order]

        print(comb_cov)
        exit()

        # Remove old dataframes.
        del eigen_df

        return comb_cov

    def save(self):
        save_dataframe(df=self.cov_df, outpath=self.outpath,
                       index=True, header=True)

    def clear_variables(self):
        self.tech_covs_file = None
        self.tech_covs_df = None
        self.cohort_file = None
        self.cohort_df = None
        self.eig_file = None
        self.n_eigen = None
        self.decon_file = None
        self.decon_df = None
        self.sample_dict = None
        self.sample_order = None
        self.force = None

    def print_arguments(self):
        print("Arguments:")
        if self.tech_covs_df is not None:
            print("  > Technival covariates input shape: {}".format(self.tech_covs_df.shape))
        else:
            print("  > Technical covariates input file: {}".format(self.tech_covs_file))
        if self.cohort_df is not None:
            print("  > Cohort input shape: {}".format(self.cohort_df.shape))
        else:
            print("  > Cohort input file: {}".format(self.cohort_file))
        print("  > Eigenvectors input file: {}".format(self.eig_file))
        print("  > N. Eigenvectors: {}".format(self.n_eigen))
        if self.decon_df is not None:
            print("  > Deconvolution input shape: {}".format(self.decon_df.shape))
        else:
            print("  > Deconvolution input file: {}".format(self.decon_file))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
