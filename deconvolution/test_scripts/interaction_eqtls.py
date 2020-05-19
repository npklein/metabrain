#!/usr/bin/env python3

"""
File:         interaction_eqtls.py
Created:      2020/05/01
Last Changed: 2020/05/19
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
import pandas as pd
import matplotlib

matplotlib.use('Agg')

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
        # Interaction matrix files.
        category = 'ms_output'
        self.inter_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/interaction_table.txt.gz".format(
            category)
        self.tvalue_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/{}/covariates/tvalue_table.txt.gz".format(
            category)
        self.eqtl_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/combine_eqtlprobes/eQTLprobes_combined.txt.gz".format(
            category)
        self.geno_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/{}/create_matrices/genotype_table.txt.gz".format(
            category)
        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'interaction_eqtls', category)
        self.cols_of_interest = ["SEX", "CellMapNNLS_Astrocyte", "CellMapNNLS_EndothelialCell", "CellMapNNLS_Macrophage", "CellMapNNLS_Neuron", "CellMapNNLS_Oligodendrocyte"]
        self.max_url_len = 8190

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        print("Loading dataframes.")
        inter_df = pd.read_csv(self.inter_path, sep="\t", header=0, index_col=0)
        print("\tInteraction matrix: {}".format(inter_df.shape))

        tvalue_df = pd.read_csv(self.tvalue_path, sep="\t", header=0, index_col=0)
        print("\tT-value matrix: {}".format(tvalue_df.shape))

        eqtl_df = pd.read_csv(self.eqtl_path, sep="\t", header=0)
        print("\teQTL matrix: {}".format(eqtl_df.shape))

        geno_df = pd.read_csv(self.geno_path, sep="\t", header=0, index_col=0)
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
        print("Determining flip column.")
        flip_col = []
        for i in range(len(eqtl_df.index)):
            flip = 1
            genotype = geno_df.iloc[i, :].round(0)
            counts = genotype.value_counts()
            if -1 in counts.index:
                counts = counts.drop([-1], axis=0)
            if counts.idxmin() != 2.0:
                flip = -1
            flip_col.append(flip)

        print(pd.Series(flip_col).value_counts())

        print("Working.")
        for col in self.cols_of_interest:
            print("Working on column: {}".format(col))

            pos_filename = os.path.join(self.outdir,
                                        "{}_positive.txt".format(col))
            neg_filename = os.path.join(self.outdir,
                                        "{}_negative.txt".format(col))

            df = eqtl_df.copy()
            df["zscore"] = inter_df.loc[col, :].values
            df["tvalue"] = (tvalue_df.loc[col, :].values * flip_col)

            df = df.loc[df.loc[:, "zscore"] > 2.49, :]

            positive = df.loc[df.loc[:, "tvalue"] > 0, :].copy()
            positive.sort_values(by="zscore", ascending=False, inplace=True)
            self.save(positive, pos_filename, self.max_url_len)

            negative = df.loc[df.loc[:, "tvalue"] < 0, :].copy()
            negative.sort_values(by="zscore", ascending=False, inplace=True)
            self.save(negative, neg_filename, self.max_url_len)

    @staticmethod
    def save(df, outfile, max_url_len):
        if len(df.index) <= 0:
            return

        print("\tWriting output file: {}".format(outfile))
        with open(outfile, 'w') as f:
            f.write("Index\tSNPName\tProbeName\tHGNCName\tZ-score\tT-value\n")

            url_string = ""
            url_genes = []
            for i, row in df.iterrows():
                f.write("{}\t{}\t{}\t{}\t"
                        "{:.2f}\t{:.2f}\n".format(i,
                                                  row["SNPName"],
                                                  row["ProbeName"],
                                                  row["HGNCName"],
                                                  row["zscore"],
                                                  row["tvalue"],
                                                  ))
                if (len(url_string) + len(row["ProbeName"])) < max_url_len:
                    if row["HGNCName"] not in url_genes:
                        url_string += row["ProbeName"]
                        url_genes.append(row["HGNCName"])

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
