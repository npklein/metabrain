#!/usr/bin/env python3

"""
File:         interaction_eqtls.py
Created:      2020/05/01
Last Changed: 2020/05/20
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
import os

# Third party imports.
import scipy.stats as st
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Interaction eQTLs"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class main():
    def __init__(self):
        self.input = {'cis_output': 6.31e-03,
                      'alzheimer_output': 1.35e-03,
                      'schizophrenia_output': 1.63e-03,
                      'parkinson_output': 2.24e-03,
                      'depression_output': 1.52e-03,
                      'ms_output': 5.71e-03,
                      'alzheimer_trans_output': 1.09e-03,
                      'schizophrenia_trans_output': 1.24e-02,
                      'parkinson_trans_output': 3.10e-03,
                      'depression_trans_output': 1.86e-03,
                      'ms_trans_output': 8.62e-04}
        self.cols_of_interest = ["SEX", "CellMapNNLS_Astrocyte",
                                 "CellMapNNLS_EndothelialCell",
                                 "CellMapNNLS_Macrophage", "CellMapNNLS_Neuron",
                                 "CellMapNNLS_Oligodendrocyte"]
        self.maf_cutoff = 0.05
        self.max_url_len = 8190
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'interaction_eqtls')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        for dir, pval_cutoff in self.input.items():
            z_score_cutoff = abs(st.norm.isf(pval_cutoff))

            print("Directory: {}".format(dir))
            print("Z-score cutoff: {:.2f}".format(z_score_cutoff))

            # Define filenames.
            inter_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/interaction_table.txt.gz".format(
                dir)
            tvalue_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/tvalue_table.txt.gz".format(
                dir)
            eqtl_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/combine_eqtlprobes/eQTLprobes_combined.txt.gz".format(
                dir)
            geno_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/create_matrices/genotype_table.txt.gz".format(
                dir)
            outdir = os.path.join(self.outdir, dir)

            if not os.path.exists(outdir):
                os.makedirs(outdir)

            print("Loading dataframes.")
            inter_df = pd.read_csv(inter_path, sep="\t", header=0, index_col=0)
            print("\tInteraction matrix: {}".format(inter_df.shape))

            tvalue_df = pd.read_csv(tvalue_path, sep="\t", header=0,
                                    index_col=0)
            print("\tT-value matrix: {}".format(tvalue_df.shape))

            eqtl_df = pd.read_csv(eqtl_path, sep="\t", header=0)
            print("\teQTL matrix: {}".format(eqtl_df.shape))

            geno_df = pd.read_csv(geno_path, sep="\t", header=0, index_col=0)
            print("\tgenotype matrix: {}".format(geno_df.shape))

            # Check if the files math up.
            for i in range(inter_df.shape[1]):
                snp_name = eqtl_df["SNPName"][i]
                if not inter_df.columns[i].startswith(snp_name) or not \
                        tvalue_df.columns[i].startswith(snp_name):
                    print("Input files do not match (1).")
                    exit()
                if not geno_df.index[i] == snp_name:
                    print("Input files do not match (2).")
                    exit()

            # Create a flip column.
            print("Getting genotype info.")
            eqtl_geno_info = []
            for i in range(len(eqtl_df.index)):
                flip = 1
                genotype = geno_df.iloc[i, :].round(0)
                counts = genotype.value_counts()
                for x in [0.0, 1.0, 2.0]:
                    if x not in counts:
                        counts.loc[x] = 0

                zero_geno_count = (counts[0.0] * 2) + counts[1.0]
                two_geno_count = (counts[2.0] * 2) + counts[1.0]
                if two_geno_count > zero_geno_count:
                    flip = -1

                sum = zero_geno_count + two_geno_count
                minor_allele_frequency = min(zero_geno_count,
                                             two_geno_count) / sum

                eqtl_geno_info.append(
                    [flip, sum, counts[0.0], counts[1.0], counts[2.0],
                     minor_allele_frequency])
            eqtl_geno_info_df = pd.DataFrame(eqtl_geno_info,
                                             columns=["flip", "N", "0.0", "1.0",
                                                      "2.0", "MAF"])
            print("Flips:")
            print(eqtl_geno_info_df["flip"].value_counts())
            print("MAF below cutoff:")
            print(len(eqtl_geno_info_df.loc[
                      eqtl_geno_info_df["MAF"] < self.maf_cutoff, :].index))

            # Combine the data and drop the results with too few samples.
            df = eqtl_df.merge(eqtl_geno_info_df, left_index=True,
                               right_index=True)

            print("Working on column:")
            for col in self.cols_of_interest:
                print("\t{}".format(col))

                pos_filename = os.path.join(outdir,
                                            "{}_positive.txt".format(col))
                neg_filename = os.path.join(outdir,
                                            "{}_negative.txt".format(col))

                working_df = df.copy()
                working_df["zscore"] = inter_df.loc[col, :].values
                working_df["raw_tvalue"] = tvalue_df.loc[col, :].values
                working_df["tvalue"] = working_df["raw_tvalue"] * working_df[
                    "flip"]

                working_df = working_df.loc[(working_df.loc[:,
                                             "MAF"] > self.maf_cutoff) &
                                            (working_df.loc[:,
                                             "zscore"] > z_score_cutoff), :]
                if len(working_df.index) <= 0:
                    continue

                positive = working_df.loc[working_df.loc[:, "tvalue"] > 0,
                           :].copy()
                positive.sort_values(by="zscore", ascending=False, inplace=True)
                self.save(positive, pos_filename, self.max_url_len,
                          z_score_cutoff)

                negative = working_df.loc[working_df.loc[:, "tvalue"] < 0,
                           :].copy()
                negative.sort_values(by="zscore", ascending=False, inplace=True)
                self.save(negative, neg_filename, self.max_url_len,
                          z_score_cutoff)

    @staticmethod
    def save(df, outfile, max_url_len, z_score_cutoff):
        if len(df.index) <= 0:
            return

        print("\tWriting output file: {}".format(outfile))
        with open(outfile, 'w') as f:
            f.write(
                "Index\tSNPName\tProbeName\tHGNCName\tN\tMAF\tZ-score\tT-value\n")

            url_string = ""
            url_genes = []
            for i, row in df.iterrows():
                f.write("{}\t{}\t{}\t{}\t{}\t{:.2f}\t"
                        "{:.2f}\t{:.2f}\n".format(i,
                                                  row["SNPName"],
                                                  row["ProbeName"],
                                                  row["HGNCName"],
                                                  row["N"],
                                                  row["MAF"],
                                                  row["zscore"],
                                                  row["tvalue"],
                                                  ))
                if (len(url_string) + len(row["ProbeName"])) < max_url_len:
                    if row["HGNCName"] not in url_genes:
                        url_string += row["ProbeName"]
                        url_genes.append(row["HGNCName"])

            f.write("\n")

            f.write("Z-score cutoff: {}\n".format(z_score_cutoff))
            f.write("\n")

            unique_snps = set(df["SNPName"])
            f.write("Number of unique SNPs: {}\n".format(len(unique_snps)))
            f.write("Unique SNPs:\n{}\n".format(', '.join(unique_snps)))
            f.write("\n")

            unique_genes = set(df["HGNCName"])
            f.write("Number of unique genes: {}\n".format(len(unique_genes)))
            f.write("Unique genes:\n{}\n".format(', '.join(unique_genes)))
            f.write("\n")

            f.write("Number of URL genes: {}\n".format(len(url_genes)))
            f.write("URL genes:\n{}\n".format(', '.join(url_genes)))
            f.write("\n")

        f.close()


if __name__ == '__main__':
    m = main()
    m.start()
