#!/usr/bin/env python3

"""
File:         plot_missingness.py
Created:      2020/03/04
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
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.


# Metadata.
__program__ = "Plot Missingness"
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

    def __init__(self, geno_file, cov_file, cohorts):
        """
        Initializer method for the main class.

        :param geno_file: string, the genotype data input file.
        :param cov_file: string, the covariates data input file.
        :param cohorts: list, the cohorts columns.
        """
        self.geno_file = geno_file
        self.cov_file = cov_file
        self.cohorts = cohorts
        self.outdir = os.path.join(os.getcwd(), 'output', 'plots')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        sns.set()

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for developme
        """
        # Print arguments.
        print("Arguments:")
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Covariate input file: {}".format(self.cov_file))
        print("")

        # Load the genotype data.
        print("\tLoading genotype matrix.")
        geno_df = pd.read_csv(self.geno_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        n_snps = geno_df.shape[0]
        print("\t\tShape: {}".format(geno_df.shape))

        # Load the coviarate data.
        print("\tLoading covariate matrix.")
        cov_df = pd.read_csv(self.cov_file, sep="\t", header=0, index_col=0,
                             nrows=nrows)
        print("\t\tShape: {}".format(cov_df.shape))


        # Subset the cohort.
        counts = pd.DataFrame(0, index=self.cohorts, columns=["samples", "snps", "complete", "incomplete"])
        for cohort in self.cohorts:
            counts.at[cohort, "snps"] = n_snps

            # Get the samples of that cohort.
            samples = cov_df.loc[cohort, :].to_frame()
            samples = samples.loc[samples[cohort] == 1, :]
            n_samples = samples.shape[0]
            counts.at[cohort, "samples"] = n_samples

            # Subset the genotype dataframe.
            subset = geno_df.loc[:, samples.index]

            # Count number of -1.
            subset.replace(-1, np.nan, inplace=True)
            mask = ~subset.isnull().any(axis=1)
            mask.reset_index(drop=True, inplace=True)
            n_complete = mask.value_counts()[True]
            n_incomplete = subset.shape[0] - n_complete

            counts.at[cohort, "complete"] = n_complete
            counts.at[cohort, "incomplete"] = n_incomplete
            counts.at[cohort, "prcnt"] = round((100 / n_snps) * n_complete, 2)

        print("Results")
        print(counts)


if __name__ == "__main__":
    GENOTYPE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                            "output", "2019-11-06-FreezeTwoDotOne",
                            "2020-03-03-interaction-analyser",
                            "step4-prepare-matrices", "output",
                            "unmasked", "genotype_table.txt.gz")

    COVARIATES = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output", "2019-11-06-FreezeTwoDotOne",
                              "2020-03-03-interaction-analyser",
                              "step4-prepare-matrices", "output",
                              "unmasked", "covariate_table.txt.gz")

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
               "UCLA_ASD-EUR",
               "ENA-EU"]

    # Start the program.
    MAIN = Main(geno_file=GENOTYPE,
                cov_file=COVARIATES,
                cohorts=COHORTS)

    MAIN.start()
