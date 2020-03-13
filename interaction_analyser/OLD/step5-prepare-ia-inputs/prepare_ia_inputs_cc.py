#!/usr/bin/env python3

"""
File:         prepare_ia_inputs_cc.py
Created:      2020/03/04
Last Changed: 2020/03/05
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


# Metadata.
__program__ = "Prepare Interaction Analyser Inputs (Complete Case)"
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

    def __init__(self, geno_file, expr_file, cov_file, exclude, eqtl_ia):
        """
        Initializer method for the main class.

        :param geno_file: string, the genotype data input file.
        :param expr_file: string, the expression data input file.
        :param cov_file: string, the covariates data input file.
        :param exclude: list, the cohorts to exclude.
        :param eqtl_ia: string, the path to the interaction analyser.
        """
        self.geno_file = geno_file
        self.expr_file = expr_file
        self.cov_file = cov_file
        self.exclude = exclude
        self.eqtl_ia = eqtl_ia

        self.cc_dir = os.path.join(os.getcwd(), 'output_cc', 'complete_case')
        self.binary_dir = os.path.join(os.getcwd(), 'output_cc', 'binary')

        for dir in [self.cc_dir, self.binary_dir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > Covariate input file: {}".format(self.cov_file))
        print("  > Exclude cohorts: {}".format(self.exclude))
        print("  > eQTL Interaction Analyser: {}".format(self.eqtl_ia))
        print("  > Complete-case output directory: {}".format(self.cc_dir))
        print("  > Binary output directory: {}".format(self.binary_dir))
        print("")

        # Load the genotype data.
        print("\tLoading genotype matrix.")
        geno_df = pd.read_csv(self.geno_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        print("\t\tShape: {}".format(geno_df.shape))

        # Load the expression data.
        print("\tLoading expression matrix.")
        expr_df = pd.read_csv(self.expr_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        print("\t\tShape: {}".format(expr_df.shape))

        # Load the coviarate data.
        print("\tLoading covariate matrix.")
        cov_df = pd.read_csv(self.cov_file, sep="\t", header=0, index_col=0,
                             nrows=nrows)
        print("\t\tShape: {}".format(cov_df.shape))

        # Create a mask on which to filter the SNPs.
        mask = self.create_mask(cov_df, geno_df, self.exclude)

        print("Working on covariate matrix")
        cov_cc_path = os.path.join(self.cc_dir,
                                   os.path.basename(self.cov_file).split(".")[
                                       0] + '_cc.txt')
        self.write_uncompressed_file(cov_df, cov_cc_path)
        self.create_binary_file(cov_cc_path, 'Covariates')
        del cov_df

        # Start with the genotype.
        print("Working on genotype matrix")
        print("\tShape before: {}".format(geno_df.shape))
        geno_df = geno_df.loc[mask, :]
        print("\tShape after: {}".format(geno_df.shape))
        geno_cc_path = os.path.join(self.cc_dir,
                                    os.path.basename(self.geno_file).split(".")[
                                        0] + '_cc.txt')
        self.write_uncompressed_file(geno_df, geno_cc_path)
        self.create_binary_file(geno_cc_path, 'Genotypes')
        del geno_df

        print("Working on expression matrix")
        print("\tShape before: {}".format(expr_df.shape))
        expr_df = expr_df.loc[mask, :]
        print("\tShape after: {}".format(expr_df.shape))
        expr_cc_path = os.path.join(self.cc_dir,
                                    os.path.basename(self.expr_file).split(".")[
                                        0] + '_cc.txt')
        self.write_uncompressed_file(expr_df, expr_cc_path)
        self.create_binary_file(expr_cc_path, 'Expression')
        del expr_df

    @staticmethod
    def create_mask(cov_df, geno_df, exclude):
        """
        Method for creating a mask for the dataframes.

        :param df: data frame, the input data frame.
        :return mask: boolean array, the mask of missing genotypes.nt.
        """
        # Create mask for the samples we are interested in.
        sample_mask = pd.Series(True, index=geno_df.columns)
        if exclude is not None:
            for cohort in exclude:
                # Get the samples of that cohort.
                samples = cov_df.loc[cohort, :].copy()
                samples = samples.to_frame()
                sub_mask = samples[cohort] == 0
                sample_mask = (sample_mask & sub_mask)

        sample_counts = sample_mask.value_counts()
        including_samples = sample_counts[True]
        excluding_samples = cov_df.shape[1] - including_samples
        print("Including {} samples \tExcluding {} "
              "samples.".format(including_samples, excluding_samples))

        # Subset those samples in the genotype matrix.
        subset = geno_df.loc[:, sample_mask].copy()
        subset.replace(-1, np.nan, inplace=True)
        mask = ~subset.isnull().any(axis=1)

        snp_counts = mask.value_counts()
        including_snps = snp_counts[True]
        excluding_snps = subset.shape[0] - including_snps
        print("Including {} SNP's \tExcluding {} "
              "SNP's.".format(including_snps, excluding_snps))

        return mask

    @staticmethod
    def write_uncompressed_file(df, outpath):
        """
        Method for writing a uncompressed matrix file.

        :param df: dataframe, the input pandas dataframe.
        :param outpath: string, the path for the output file.
        :return: outpath, string, the output file path.
        """
        # Write output file.
        print("\tWriting uncompressed file.")
        df.to_csv(outpath, sep="\t")

        return outpath

    def create_binary_file(self, inpath, out_filename):
        """
        Method for creating a binary file for the Interaction Analyser input.

        :param inpath: string, the input file path.
        :param out_filename: string, the output filename.
        :return:
        """
        print("\tCreating binary file.")
        outpath = os.path.join(self.binary_dir, out_filename + '.binary')
        print(outpath)

        # Convert to binary.
        command = 'java -jar {} --convertMatrix -i {} -o {}'.format(
            self.eqtl_ia, inpath, outpath)
        print("\t\t{}".format(command))
        os.system(command)


if __name__ == "__main__":
    GENOTYPE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                            "output", "2019-11-06-FreezeTwoDotOne",
                            "2020-03-03-interaction-analyser",
                            "step4-prepare-matrices", "output",
                            "unmasked", "genotype_table.txt.gz")

    EXPRESSION = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output", "2019-11-06-FreezeTwoDotOne",
                              "2020-03-03-interaction-analyser",
                              "step4-prepare-matrices", "output",
                              "unmasked", "expression_table.txt.gz")

    COVARIATES = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output", "2019-11-06-FreezeTwoDotOne",
                              "2020-03-03-interaction-analyser",
                              "step4-prepare-matrices", "output",
                              "unmasked", "covariate_table.txt.gz")

    EQTL_IA = os.path.join(os.path.sep, "groups", "umcg-biogen",
                           "tmp03", "output", "2019-11-06-FreezeTwoDotOne",
                           "2020-03-03-interaction-analyser",
                           "step6-interaction-analyser",
                           "eQTLInteractionAnalyser-1.2-SNAPSHOT-jar-with-dependencies.jar")

    #EXCLUDE = ["ENA-EU"]
    EXCLUDE = None

    # Start the program.
    MAIN = Main(geno_file=GENOTYPE,
                expr_file=EXPRESSION,
                cov_file=COVARIATES,
                exclude=EXCLUDE,
                eqtl_ia=EQTL_IA)

    MAIN.start()
