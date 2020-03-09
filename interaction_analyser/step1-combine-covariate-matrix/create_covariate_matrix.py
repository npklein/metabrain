#!/usr/bin/env python3

"""
File:         create_covariate_matrix.py
Created:      2020/02/25
Last Changed: 2020/03/09
Author:       M.Vochteloo

Copyright (C) 2019 M.Vochteloo

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


# Metadata.
__program__ = "Create Covariate Matrix"
__author__ = "M. Vochteloo"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class Main:
    """
    Main class of the program.
    """

    def __init__(self, cov_file, tech_covs, cohorts, pheno_file, eig_file,
                 n_comps, out_filename):
        """
        The initializer for the main class.

        :param cov_file: string, the covariate data input file.
        :param tech_covs: list, the technical coviarates columns
        :param cohorts: list, the cohorts columns.
        :param pheno_file: string, the phenotype data input file.
        :param eig_file: string, the eigenvectors data input file.
        :param n_comps: int, the number of components to include.
        :param out_filename: string, the output filename.
        """
        self.cov_file = cov_file
        self.tech_covs = tech_covs
        self.cohorts = cohorts
        self.pheno_file = pheno_file
        self.eig_file = eig_file
        self.n_comps = n_comps
        outdir = os.path.join(os.getcwd(), 'output')
        self.outpath = os.path.join(outdir, out_filename)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def start(self):
        """
        Main method for the main class. Does all the work.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Covariate input file: {}".format(self.cov_file))
        print("  > Technical covariate columns: {}".format(self.tech_covs))
        print("  > Cohort columns: {}".format(self.cohorts))
        print("  > Phenotype input file: {}".format(self.pheno_file))
        print("  > Eigenvectors input file: {}".format(self.eig_file))
        print("  > N components to include: {}".format(self.n_comps))
        print("  > Output path: {}".format(self.outpath))
        print("")

        # read the covariates file.
        print("Loading covariate matrix.")
        cov_df = pd.read_csv(self.cov_file, sep="\t", header=0, index_col=0)
        print("\tShape: {}".format(cov_df.shape))

        # filter thechnical covariates and the cohorts out.
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
        pheno_df = pd.read_csv(self.pheno_file, sep="\t", header=0, index_col=4,
                               low_memory=False)
        pheno_df = pheno_df.loc[:, ["Gender", "sex.by.expression"]]
        pheno_df.replace("no expression available", np.nan, inplace=True)
        pheno_df["GENDER"] = pheno_df['sex.by.expression'].combine_first(pheno_df['Gender'])
        gender_df = pheno_df["GENDER"].to_frame()
        del pheno_df
        gender_df = gender_df.replace({"GENDER": {"M": 0, "F": 1, np.nan: -1}})
        print("\tShape: {}".format(gender_df.shape))

        # read the eigenvectors file.
        print("Loading eigenvectors matrix.")
        eigen_df = pd.read_csv(self.eig_file, sep="\t", header=0, index_col=0)
        eigen_df = eigen_df.loc[:,
                   ["Comp{}".format(x) for x in range(1, self.n_comps + 1)]]
        print("\tShape: {}".format(eigen_df.shape))

        # merge.
        print("Merging matrices.")
        tmp_df1 = tech_cov_df.merge(cohorts_df, left_index=True,
                                    right_index=True)
        tmp_df2 = tmp_df1.merge(gender_df, left_index=True, right_index=True)
        merged_df = tmp_df2.merge(eigen_df, left_index=True, right_index=True)
        del tmp_df1, tmp_df2
        # print(merged_df.iloc[0, :])
        print("\tShape: {}".format(merged_df.shape))

        with pd.option_context('display.max_rows', None, 'display.max_columns',
                               None):  # more options can be specified also
            print(merged_df.loc["AN11864_ba41.42.22"])

        # post process the matrix.
        print("Post-process covariate matrix.")
        merged_df = merged_df.T
        merged_df.index = merged_df.index.set_names(['Sample'])
        print("\tShape: {}".format(merged_df.shape))

        # write outfile.
        merged_df.to_csv(self.outpath, sep="\t", compression='gzip')


if __name__ == "__main__":
    # Define main variables.
    COVARIATES = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output",
                              "2019-11-06-FreezeTwoDotOne",
                              "2020-01-31-expression-tables",
                              "2020-02-05-step6-covariate-removal",
                              "2020-02-17-step5-remove-covariates-per-dataset",
                              "2020-02-17-freeze2dot1.TMM.Covariates.withBrainRegion-noncategorical-variable.top20correlated-cortex-withMDS.txt.gz"
                              )

    TECHNICAL_COVS = ["PCT_MRNA_BASES",
                      "PCT_INTRONIC_BASES",
                      "MEDIAN_3PRIME_BIAS",
                      "PCT_USABLE_BASES",
                      "PCT_INTERGENIC_BASES",
                      "PCT_UTR_BASES",
                      "PCT_READS_ALIGNED_IN_PAIRS",
                      "PCT_CHIMERAS",
                      "PF_READS_IMPROPER_PAIRS",
                      "PF_HQ_ALIGNED_Q20_BASES",
                      "PF_HQ_ALIGNED_BASES",
                      "PCT_PF_READS_IMPROPER_PAIRS",
                      "PF_READS_ALIGNED",
                      "avg_mapped_read_length",
                      "avg_input_read_length",
                      "uniquely_mapped",
                      "total_reads",
                      "Total.Sequences_R1",
                      "MDS1",
                      "MDS2",
                      "MDS3",
                      "MDS4"]

    COHORTS = ["AMPAD-MSBB-V2-AFR",
               "CMC-AFR",
               "LIBD_1M-AFR",
               "LIBD_h650-AFR",
               "AMPAD-MAYO-V2-EUR",
               "AMPAD-MSBB-V2-EUR",
               "AMPAD-ROSMAP-V2-EUR",
               "BrainGVEX-V2-EUR",
               "CMC-EUR",
               "GTEx-EUR",
               "GVEx",
               "LIBD_1M-EUR",
               "LIBD_h650-EUR",
               "NABEC-H550-EUR",
               "NABEC-H610-EUR",
               "TargetALS-EUR",
               "UCLA_ASD-EUR"]

    PHENOTYPE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                             "output",
                             "2019-11-06-FreezeTwoDotOne",
                             "2020-02-03-phenotype-table",
                             "2020-03-03.brain.phenotypes.txt")

    GENDER_COLUMN = "sex.by.expression"

    EIGENVECTORS = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                                "output",
                                "2019-11-06-FreezeTwoDotOne",
                                "2020-01-31-expression-tables",
                                "2020-02-05-step6-covariate-removal",
                                "2020-02-17-step5-remove-covariates-per-dataset",
                                "output-cortex",
                                "MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.PCAOverSamplesEigenvectors.txt.gz")

    NUM_OF_COMPONENTS = 50

    OUT_FILENAME = "covariates-cortex-withMD-withGender-with{}PCs.txt.gz".format(
        NUM_OF_COMPONENTS)

    # Start the program.
    MAIN = Main(cov_file=COVARIATES,
                tech_covs=TECHNICAL_COVS,
                cohorts=COHORTS,
                pheno_file=PHENOTYPE,
                eig_file=EIGENVECTORS,
                n_comps=NUM_OF_COMPONENTS,
                out_filename=OUT_FILENAME)
    MAIN.start()
