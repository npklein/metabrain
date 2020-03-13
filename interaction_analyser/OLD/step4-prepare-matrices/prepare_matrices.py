#!/usr/bin/env python3

"""
File:         prepare_matrices.py
Created:      2020/02/25
Last Changed: 2020/03/03
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
import argparse
import gzip
import os

# Third party imports.
import pandas as pd

# Local application imports.


# Metadata.
__program__ = "Prepare Matrices"
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

    def __init__(self, gte_file, eqtl_file, cov_file, geno_file, expr_file,
                 masked):
        """
        Initializer method for the main class.

        :param gte_file: string, the translation file for the samples.
        :param eqtl_file: string, the file containing the eQTL effects of
                          interest + also acts as translation file.
        :param cov_file: string, the covariates data input file.
        :param geno_file: string, the genotype data input file.
        :param expr_file: string, the expression data input file.
        :param masked: boolean, whether or not to mask the rows / columns of
                       the output matrices.
        """
        self.gte_file = gte_file
        self.eqtl_file = eqtl_file
        self.cov_file = cov_file
        self.geno_file = geno_file
        self.expr_file = expr_file
        self.masked = masked

        if self.masked:
            self.outdir = os.path.join(os.getcwd(), 'output', 'masked')
        else:
            self.outdir = os.path.join(os.getcwd(), 'output', 'unmasked')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Print arguments.
        print("Arguments:")
        print("  > GTE input  file: {}".format(self.gte_file))
        print("  > eQTL input file: {}".format(self.eqtl_file))
        print("  > Covariate input file: {}".format(self.cov_file))
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > Masking SNP's / samples: {}".format(self.masked))
        print("  > Output directory: {}".format(self.outdir))
        print("")

        # Define the names of the output matrices.
        geno_outfile = os.path.join(self.outdir, "genotype_table.txt.gz")
        expr_outfile = os.path.join(self.outdir, "expression_table.txt.gz")
        cov_outfile = os.path.join(self.outdir, "covariate_table.txt.gz")
        allele_outfile = os.path.join(self.outdir, "genotype_alleles.txt.gz")

        # Remove the output files if they already exist.
        print("Preparing output files.")
        for outfile in [geno_outfile, expr_outfile, cov_outfile,
                        allele_outfile]:
            if os.path.isfile(outfile):
                os.remove(outfile)

        # Create the sample translation data frame.
        print("Creating sample translation data frame.")
        gte_df = pd.read_csv(self.gte_file, sep="\t", header=None)
        if self.masked:
            gte_df["masked"] = ["sample_" + str(x) for x in
                                range(gte_df.shape[0])]
            gte_df.to_csv(
                os.path.join(self.outdir, "sample_translate_table.txt.gz"),
                sep="\t", index=False, compression='gzip')
        sample_dict, sample_mask = self.create_sample_map(gte_df)
        del gte_df

        # Get the SNP translation data frame.
        print("Loading eQTL matrix.")
        eqtl_df = pd.read_csv(self.eqtl_file, sep="\t", header=0)
        if self.masked:
            subset_eqtl = eqtl_df.loc[:, ["SNPName", "ProbeName"]]
            subset_eqtl["masked"] = ["SNP_" + str(x) for x in
                                     range(subset_eqtl.shape[0])]
            subset_eqtl.to_csv(os.path.join(self.outdir,
                                            "SNP_translate_table.txt.gz"),
                               sep="\t", index=False, compression='gzip')
        n_snps = eqtl_df.shape[0]
        print("\tShape: {}".format(eqtl_df.shape))

        # Load the covariate data frame.
        print("Loading covariate matrix.")
        cov_df = pd.read_csv(self.cov_file, sep="\t", header=0, index_col=0)
        print("\tShape: {}".format(cov_df.shape))

        # Construct the covariate translation data frame.
        print("Creating covariate translation data frame.")
        cov_trans_df = pd.DataFrame(cov_df.index)
        if self.masked:
            cov_trans_df["masked"] = ["cov_" + str(x) for x in
                                      range(cov_trans_df.shape[0])]
            cov_trans_df.to_csv(
                os.path.join(self.outdir, "covariate_translate_table.txt.gz"),
                sep="\t", index=False, compression='gzip')
        cov_dict = self.create_cov_dict(cov_trans_df)

        # Rename and reorder the covariate matrix.
        cov_df = cov_df.rename(index=cov_dict, columns=sample_dict)
        cov_df = cov_df.loc[:, sample_mask]
        cov_df.to_csv(cov_outfile, sep="\t", compression='gzip')
        del cov_df, cov_dict

        # Load the genotype matrix file.
        print("Loading genotype matrix.")
        geno_df = pd.read_csv(self.geno_file, sep="\t", header=0, nrows=nrows)
        allele_df = geno_df.loc[:, ["SNP", "Alleles", "MinorAllele"]].copy()
        geno_df = geno_df.drop(["Alleles", "MinorAllele"], axis=1)
        print("\tShape: {}".format(geno_df.shape))

        # Load the expression matrix file.
        print("Loading expression matrix.")
        expr_df = pd.read_csv(self.expr_file, sep="\t", header=0, nrows=nrows)
        print("\tShape: {}".format(expr_df.shape))

        # Construct the genotype / expression matrices.
        print("Constructing matrices.")
        geno_str_buffer = ["-" + "\t" + "\t".join(sample_mask) + "\n"]
        expr_str_buffer = ["-" + "\t" + "\t".join(sample_mask) + "\n"]
        allele_str_buffer = ["\t".join(list(allele_df.columns)) + "\n"]
        for i, row in eqtl_df.iterrows():
            if (i % 250 == 0) or (i == (n_snps - 1)):
                print("\t Processing {}/{} [{:.2f}%]".format(i, n_snps,
                                                             (100 / n_snps) * i))

                # Write output files.
                self.write_buffer(geno_outfile, geno_str_buffer)
                geno_str_buffer = []

                self.write_buffer(expr_outfile, expr_str_buffer)
                expr_str_buffer = []

                self.write_buffer(allele_outfile, allele_str_buffer)
                allele_str_buffer = []

            # Get the eQTL
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            index = "SNP_" + str(i)
            if not self.masked:
                index = snp_name

            # Used for development.
            # snp_name = "10:100145864:rs4919426:T_C"
            # probe_name = "ENSG00000000003.15"
            # End used for development.

            # Get the genotype.
            genotype = geno_df.loc[geno_df["SNP"] == snp_name, :]
            genotype = genotype.rename(columns=sample_dict)
            genotype = genotype.loc[:, sample_mask]
            if (len(genotype.index)) != 1:
                print("SNP: {} gives 0 or >1 "
                      "genotypes.".format(snp_name))
                continue
            geno_str = index + "\t" + "\t".join(
                genotype.iloc[0, :].astype(str).values) + "\n"
            geno_str_buffer.append(geno_str)

            # Get the alleles.
            alleles = allele_df.loc[allele_df["SNP"] == snp_name, :]
            if (len(alleles.index)) != 1:
                print("SNP: {} gives 0 or >1 "
                      "alleles.".format(snp_name))
                continue
            allele_str = "{}\t{}\t{}\n".format(index,
                                               alleles.iloc[0]["Alleles"],
                                               alleles.iloc[0][
                                                   "MinorAllele"])
            allele_str_buffer.append(allele_str)

            # Get the expression.
            expression = expr_df.loc[expr_df["-"] == probe_name, :]
            expression = expression.rename(columns=sample_dict)
            expression = expression.loc[:, sample_mask]
            if (len(expression.index)) != 1:
                print("Probe: {} gives 0 or >1 expression "
                      "profiles.".format(probe_name))
                continue
            expr_str = index + "\t" + "\t".join(
                expression.iloc[0, :].astype(str).values) + "\n"
            expr_str_buffer.append(expr_str)

        # Write output files.
        if geno_str_buffer:
            self.write_buffer(geno_outfile, geno_str_buffer)

        if expr_str_buffer:
            self.write_buffer(expr_outfile, expr_str_buffer)

        if allele_str_buffer:
            self.write_buffer(allele_outfile, allele_str_buffer)

    def create_sample_map(self, df):
        """
        Method for creating a translation dict for the samples.

        :param df: pandas dataframe, the conversion dict.
        :return: dict, the translation dict.
        """
        map_dict = {}
        for _, row in df.iterrows():
            values = row.values
            keys = [values[0]]
            if self.masked:
                keys = values[:-1]
            value = values[-1]

            for key in keys:
                if key not in map_dict.keys():
                    map_dict[key] = value

        sample_order = list(df.iloc[:, -1])

        return map_dict, sample_order

    def create_cov_dict(self, df):
        """
        Method for creating a translation dict for the coviriates.

        :param df: pandas dataframe, the conversion dict.
        :return: dict, the translation dict.
        """
        map_dict = {}
        for _, row in df.iterrows():
            values = row.values
            key = values[0]
            value = values[0]
            if self.masked:
                value = values[1]

            if key not in map_dict.keys():
                map_dict[key] = value

        return map_dict

    @staticmethod
    def write_buffer(filename, buffer):
        """
        Method for writing a list of strings to a gzipped file,

        :param filename: string, the output file path.
        :param buffer: list, the lines of strings to write.
        """
        # Write output files.
        if os.path.isfile(filename):
            mode = 'ab'
        else:
            mode = 'wb'

        with gzip.open(filename, mode) as f:
            for line in buffer:
                f.write(line.encode())
        f.close()


class CommandLineArguments:
    """
    Class for the command-line arguments.
    """

    def __init__(self, program, description, version):
        """
        Initializer for the command line arguments class.

        :param program: string, the name of the program.
        :param description: string, the description of the program.
        :param version: float, the version of the program.
        """
        # Safe variables.
        self.program = program
        self.description = description
        self.version = version

        # Get the arguments.
        parser = self.create_argument_parser()
        self.arguments = parser.parse_args()

    def create_argument_parser(self):
        """
        Method for creating the ArgumentParser object.

        :return: ArgumentParser
        """
        parser = argparse.ArgumentParser(prog=self.program,
                                         description=self.description,
                                         add_help=False)
        parser.add_argument("-m",
                            "--masked",
                            action='store_true',
                            help="Whether or not to mask the rows / columns.")

        return parser

    def get_argument(self, arg_key):
        """
        Method to return the value of a command line key.

        :param arg_key: str, key in parse command line
        :return: int/str, value of arg_key. If the key did not exist,
                 return None.
        """
        if self.arguments is not None and arg_key in self.arguments:
            value = getattr(self.arguments, arg_key)
        else:
            value = None

        return value


if __name__ == "__main__":
    # Step 0 data.
    GTE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                       "output", "2019-11-06-FreezeTwoDotOne",
                       "2020-03-03-interaction-analyser",
                       "step0-combine-GTE-files", "output",
                       "GTE_combined.txt.gz")

    # Step 1 data.
    EQTLS = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-03-03-interaction-analyser",
                         "step1-combine-eQTLprobe-files", "output",
                         "eQTLProbesFDR0.05-ProbeLevel_combined.txt.gz")

    # Step 3 data.
    COVARIATE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                             "output", "2019-11-06-FreezeTwoDotOne",
                             "2020-03-03-interaction-analyser",
                             "step3-combine-covariate-matrix", "output",
                             "covariates-cortex.txt.gz")

    # Other data sources.
    GENOTYPE = os.path.join(os.path.sep, "groups", "umcg-biogen",
                            "tmp03", "output",
                            "2019-11-06-FreezeTwoDotOne",
                            "2020-02-18-eqtls",
                            "cortex-genotypedump-EURandAFR",
                            "GenotypeData.txt.gz",
                            )

    EXPRESSION = os.path.join(os.path.sep, "groups", "umcg-biogen",
                              "tmp03", "output",
                              "2019-11-06-FreezeTwoDotOne",
                              "2020-01-31-expression-tables",
                              "2020-02-05-step6-covariate-removal",
                              "2020-02-17-step5-remove-covariates-per-dataset",
                              "output-cortex",
                              "MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz"
                              )

    # Get the command line argument.
    CMA = CommandLineArguments(__program__,
                               __description__,
                               __version__)
    MASKED = CMA.get_argument("masked")

    # Start the program.
    MAIN = Main(gte_file=GTE,
                eqtl_file=EQTLS,
                cov_file=COVARIATE,
                geno_file=GENOTYPE,
                expr_file=EXPRESSION,
                masked=MASKED)

    MAIN.start()
