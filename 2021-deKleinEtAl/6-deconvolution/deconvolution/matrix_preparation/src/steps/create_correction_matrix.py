"""
File:         create_correction_matrix.py
Created:      2020/10/15
Last Changed: 2021/07/08
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
import pandas as pd
import sympy

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, \
    save_dataframe


class CreateCorrectionMatrix:
    def __init__(self, settings, log, cohort_file, cohort_df, sample_dict,
                 sample_order, force, outdir):
        self.cov_file = settings["covariates_datafile"]
        self.tech_covs = settings["technical_covariates"]
        self.mds_covs = settings["MDS_covariates"]
        self.log = log
        self.cohort_file = cohort_file
        self.cohort_df = cohort_df
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_correction_matrix')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.outpath = os.path.join(self.outdir, "correction_table.txt.gz")

        # Create empty variable.
        self.correction_df = None

    def start(self):
        self.log.info("Starting creating technical covariate file.")
        self.print_arguments()

        # Check if output file exist.
        if not check_file_exists(self.outpath) or self.force:
            self.correction_df = self.create_tech_covs_file()
            self.save()
        else:
            self.log.info("Skipping step.")

    def create_tech_covs_file(self):
        # Load the sample info.
        self.log.info("Loading covariates matrix.")
        cov_df = load_dataframe(inpath=self.cov_file,
                                header=0,
                                index_col=0,
                                logger=self.log)

        # Filter on samples and technical covariates.
        self.log.info("Filtering on samples and technical covariates.")
        cov_df.index = [self.sample_dict[x] if x in self.sample_dict else x for
                        x in cov_df.index]
        correction_df = cov_df.loc[self.sample_order, :].copy()
        del cov_df
        self.log.info("\tNew shape: {}".format(correction_df.shape))

        # loading cohort matrix.
        self.log.info("Loading cohort matrix.")
        if self.cohort_df is None:
            self.cohort_df = load_dataframe(self.cohort_file,
                                            header=0,
                                            index_col=0,
                                            logger=self.log)

        # merge.
        self.log.info("Merging matrices.")
        correction_df = pd.merge(correction_df, self.cohort_df.T, left_index=True, right_index=True)
        correction_df = correction_df.T
        correction_df.index.name = "-"

        # Validate sample order.
        if not correction_df.columns.equals(self.sample_order):
            correction_df = correction_df[self.sample_order]

        return correction_df

    def save(self):
        save_dataframe(df=self.correction_df, outpath=self.outpath,
                       index=True, header=True, logger=self.log)
        save_dataframe(df=self.correction_df.loc[self.tech_covs, :], outpath=os.path.join(self.outdir, "technical_covariates_table.txt.gz"),
                       index=True, header=True, logger=self.log)
        save_dataframe(df=self.correction_df.loc[self.mds_covs, :], outpath=os.path.join(self.outdir, "mds_covariates_table.txt.gz"),
                       index=True, header=True, logger=self.log)
        save_dataframe(df=self.correction_df.loc[self.tech_covs + self.mds_covs, :], outpath=os.path.join(self.outdir, "technical_and_mds_covariates_table.txt.gz"),
                       index=True, header=True, logger=self.log)

    def clear_variables(self):
        self.cohort_file = None
        self.cohort_df = None
        self.cov_file = None
        self.tech_covs = None
        self.mds_covs = None
        self.force = None

    def get_correction_file(self):
        return self.outpath

    def get_correction_df(self):
        return self.correction_df

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Covariates input file: {}".format(self.cov_file))
        self.log.info("  > Technical Covarates: {}".format(self.tech_covs))
        self.log.info("  > MDS Covarates: {}".format(self.mds_covs))
        if self.cohort_df is not None:
            self.log.info("  > Cohort input shape: {}".format(self.cohort_df.shape))
        else:
            self.log.info("  > Cohort input file: {}".format(self.cohort_file))
        self.log.info("  > Output path: {}".format(self.outpath))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
