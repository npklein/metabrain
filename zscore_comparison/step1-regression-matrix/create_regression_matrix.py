#!/usr/bin/env python3

"""
File:         create_regression_matrix.py
Created:      2020/02/28
Last Changed: 2020/03/04
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
import pandas as pd
from scipy import stats
import gzip

# Local application imports.


# Metadata.
__program__ = "Create Regression Matrix"
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

    def __init__(self, eqtl_file, geno_file, expr_file, allele_file,
                 out_filename):
        """
        Initializer method for the main class.

        :param eqtl_file: string, the file containing the eQTL effects of
                          interest
        :param geno_file: string, the genotype data input file.
        :param expr_file: string, the expression data input file.
        :param allele_file: string, the genotype alleles input file.
        """
        self.eqtl_file = eqtl_file
        self.geno_file = geno_file
        self.expr_file = expr_file
        self.allele_file = allele_file
        outdir = os.path.join(os.getcwd(), 'output')
        self.outpath = os.path.join(outdir, out_filename)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if os.path.isfile(self.outpath):
            os.remove(self.outpath)

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Load the eQTL file.
        print("Loading eQTL matrix.")
        eqtl_df = pd.read_csv(self.eqtl_file, sep="\t", header=0, nrows=nrows)
        nrows = eqtl_df.shape[0]
        print("\tShape: {}".format(eqtl_df.shape))

        # Load the genotype matrix file.
        print("Loading genotype matrix.")
        geno_df = pd.read_csv(self.geno_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        print("\tShape: {}".format(geno_df.shape))

        # Load the alleles
        print("Loading alleles matrix.")
        allele_df = pd.read_csv(self.allele_file, sep="\t", header=0,
                                index_col=0,
                                nrows=nrows)
        print("\tShape: {}".format(allele_df.shape))

        # Load the expression matrix file.
        print("Loading expression matrix.")
        expr_df = pd.read_csv(self.expr_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        print("\tShape: {}".format(expr_df.shape))

        # Check if shape is identical.
        if (eqtl_df.shape[0] != geno_df.shape[0]) or \
                (eqtl_df.shape[0] != expr_df.shape[0]):
            print("Input matrices rows are not identical length.")
            return
        if geno_df.shape != expr_df.shape:
            print("Genotype and expression matrices are not identical shape.")
            return

        # Prepare string buffer.
        regr_str_buffer = [
            "snp\tprobe\talleles\tminor_allele\tallele_assessed\tflipped\tslope"
            "\tintercept\tcorr_coeff\tp_value\tstd_err\toveral_z_score"
            "\tz_score_estimate\n"]

        # Correlating.
        print("Correlating:")
        for i, row in eqtl_df.iterrows():
            if (i % 250 == 0) or (i == (nrows - 1)):
                print("\t Processing {}/{} [{:.2f}%]".format(i, nrows,
                                                             (100 / nrows) * i
                                                             )
                      )

                # Write output files.
                self.write_buffer(self.outpath, regr_str_buffer)
                regr_str_buffer = []

            # Extract the usefull information.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            overal_z_score = row["OverallZScore"]
            allele_assessed = row["AlleleAssessed"]

            # Get the data.
            genotype = geno_df.iloc[i, :].T.to_frame()
            if snp_name != genotype.columns[0]:
                print("SNPName does not match in genotype subset.")
                break
            expression = expr_df.iloc[i, :].T.to_frame()
            if snp_name != expression.columns[0]:
                print("SNPName does not match in expression subset.")
                break
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]

            # Remove missing values.
            data = data.loc[(data['genotype'] >= 0.0) &
                            (data['genotype'] <= 2.0), :]

            # Determine the alleles.
            (alleles, minor_allele) = allele_df.iloc[i, :]

            # Determine whether to flip or not.
            flipped = False
            if allele_assessed != alleles.split("/")[1]:
                data['genotype'] = 2.0 - data['genotype']
                flipped = True

            # # Naive flip method.
            # if allele_assessed != minor_allele:
            #     data['genotype'] = 2.0 - data['genotype']

            # Calculate the correlation.
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                data["genotype"],
                data["expression"])

            # Calculate the z-score estimate.
            z_score_estimate = slope / std_err

            # # Naive flip method 2.0.
            # if allele_assessed != alleles.split("/")[1]:
            #     z_score_estimate = z_score_estimate * -1

            # Add to the buffer.
            regr_str = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                .format(snp_name,
                        probe_name,
                        alleles,
                        minor_allele,
                        allele_assessed,
                        flipped,
                        slope,
                        intercept,
                        r_value,
                        p_value,
                        std_err,
                        overal_z_score,
                        z_score_estimate
                        )
            regr_str_buffer.append(regr_str)

        # Write output files.
        if regr_str_buffer:
            self.write_buffer(self.outpath, regr_str_buffer)


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


if __name__ == "__main__":
    # Define main variables.
    EQTLS = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output",
                         "2019-11-06-FreezeTwoDotOne",
                         "2020-02-18-eqtls",
                         "cortex-cis-EURandAFR-iterative",
                         "Iteration1",
                         "eQTLProbesFDR0.05-ProbeLevel.txt.gz"
                         )

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

    ALLELES = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                           "output", "2019-11-06-FreezeTwoDotOne",
                           "2020-03-03-interaction-analyser",
                           "step4-prepare-matrices", "output",
                           "unmasked", "genotype_alleles.txt.gz")

    OUT_FILENAME = "regression_table.txt.gz"

    # Start the program.
    MAIN = Main(eqtl_file=EQTLS,
                geno_file=GENOTYPE,
                expr_file=EXPRESSION,
                allele_file=ALLELES,
                out_filename=OUT_FILENAME)
    MAIN.start()
