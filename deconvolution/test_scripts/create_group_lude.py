#!/usr/bin/env python3

"""
File:         create_group_lude.py
Created:      2020/03/23
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
from __future__ import print_function
from pathlib import Path
import os

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Create Group Lude"
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
        # Translate tables.
        self.cov_translate_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/mask_matrices/cov_translate_table.txt.gz"
        self.eqtl_translate_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/mask_matrices/eqtl_translate_table.txt.gz"
        self.sample_translate_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/mask_matrices/sample_translate_table.txt.gz"

        # Masked tables.
        self.eqtl_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/mask_matrices/eqtl_table.txt.gz"
        self.alleles_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/mask_matrices/genotype_alleles.txt.gz"
        self.cov1_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/mask_matrices/covariates_table.txt.gz"

        # Files from Lude.
        self.geno_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/data/genotype_table.txt.binary.CohortCorrected.txt"
        self.expr_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/data/expression_table.txt.binary.CohortAnd50PCsCorrected.txt"
        self.cov2_inpath = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/data/expression_table.txt.binary.CovariatesOptimized.txt"

        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'output')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        print("Load translate dicts.")
        eqtl_trans = self.load_translate_table(self.eqtl_translate_inpath)
        sample_trans = self.load_translate_table(self.sample_translate_inpath)
        cov_trans = self.load_translate_table(self.cov_translate_inpath)

        print("Loading Lude files.")
        geno_df, expr_df, cov2_df = self.load_lude_files()

        print("Loading masked files.")
        eqtl_df, alleles_df, cov1_df = self.load_masked_files()

        # Merge the covariate dataframe.
        cov_df = cov1_df.T.merge(cov2_df.T, left_index=True, right_index=True).T

        # Determine which samples are overlapping.
        eqtl_mask = geno_df.index
        sample_mask = cov2_df.columns

        # Subset and safe.
        print("Saving files.")
        eqtl_df = eqtl_df.iloc[[int(x.split("_")[1]) for x in eqtl_mask], :]
        eqtl_df = eqtl_df.rename(index=eqtl_trans)
        eqtl_df.to_csv(os.path.join(self.outdir, "eqtl_table.txt.gz"),
                       sep="\t", index=False, header=True, compression="gzip")
        print("\tSaved dataframe: eqtl_table.txt.gz_table with shape: "
              "{}".format(eqtl_df.shape))

        geno_df = geno_df.loc[eqtl_mask, :]
        geno_df = geno_df.loc[:, sample_mask]
        geno_df = geno_df.rename(index=eqtl_trans, columns=sample_trans)
        geno_df.index.name = "-"
        geno_df.to_csv(os.path.join(self.outdir, "genotype_table.txt.gz"),
                       sep="\t", index=True, header=True, compression="gzip")
        print("\tSaved dataframe: genotype_table.txt.gz_table with shape: "
              "{}".format(geno_df.shape))

        alleles_df = alleles_df.loc[eqtl_mask, :]
        alleles_df = alleles_df.rename(index=eqtl_trans)
        alleles_df.index.name = "-"
        alleles_df.to_csv(os.path.join(self.outdir, "genotype_alleles.txt.gz"),
                          sep="\t", index=True, header=True, compression="gzip")
        print("\tSaved dataframe: genotype_alleles.txt.gz_table with shape: "
              "{}".format(alleles_df.shape))

        expr_df = expr_df.loc[eqtl_mask, :]
        expr_df = expr_df.loc[:, sample_mask]
        expr_df = expr_df.rename(index=eqtl_trans, columns=sample_trans)
        expr_df.index.name = "-"
        expr_df.to_csv(os.path.join(self.outdir, "expression_table.txt.gz"),
                       sep="\t", index=True, header=True, compression="gzip")
        print("\tSaved dataframe: expression_table.txt.gz_table with shape: "
              "{}".format(expr_df.shape))

        cov_df = cov_df.loc[:, sample_mask]
        cov_df = cov_df.rename(index=cov_trans, columns=sample_trans)
        cov_df.index.name = "-"
        cov_df.to_csv(os.path.join(self.outdir, "covariates_table.txt.gz"),
                       sep="\t", index=True, header=True, compression="gzip")
        print("\tSaved dataframe: covariates_table.txt.gz_table with shape: "
              "{}".format(cov_df.shape))

    @staticmethod
    def load_translate_table(inpath):
        df = pd.read_csv(inpath,
                         header=0,
                         sep='\t')
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(inpath), df.shape))

        return dict(zip(df["masked"], df["unmasked"]))

    def load_lude_files(self):
        print("  Loading genotype table.")
        geno_df = pd.read_csv(self.geno_inpath,
                              header=0,
                              index_col=0,
                              sep='\t')
        geno_df_index = []
        geno_df.dropna(inplace=True)
        for key in geno_df.index:
            number = key.split("_")[1]
            geno_df_index.append("eqtl_{}".format(number))
        geno_df.index = geno_df_index
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(self.geno_inpath), geno_df.shape))

        print("  Loading expression table.")
        expr_df = pd.read_csv(self.expr_inpath,
                              header=0,
                              index_col=0,
                              sep='\t')
        expr_df.dropna(inplace=True)
        gexpr_df_index = []
        for key in expr_df.index:
            number = key.split("_")[1]
            gexpr_df_index.append("eqtl_{}".format(number))
        expr_df.index = gexpr_df_index
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(self.expr_inpath), expr_df.shape))

        print("  Loading covariate table (optimized).")
        cov2_df = pd.read_csv(self.cov2_inpath,
                              header=0,
                              index_col=0,
                              sep='\t')
        cov2_df.dropna(inplace=True)
        cov2_df.index = ["oligodendrocytes_ASPA[opti]"]
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(self.cov2_inpath), cov2_df.shape))

        return geno_df, expr_df, cov2_df

    def load_masked_files(self):
        print("  Loading eQTL table.")
        eqtl_df = pd.read_csv(self.eqtl_inpath,
                              header=0,
                              index_col=None,
                              sep='\t')
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(self.eqtl_inpath), eqtl_df.shape))

        print("  Loading alleles table.")
        alleles_df = pd.read_csv(self.alleles_inpath,
                                 header=0,
                                 index_col=0,
                                 sep='\t')
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(self.alleles_inpath),
                          alleles_df.shape))

        print("  Loading covariate table.")
        cov1_df = pd.read_csv(self.cov1_inpath,
                              header=0,
                              index_col=0,
                              sep='\t')
        print("\tLoaded dataframe: {} with shape: "
              "{}".format(os.path.basename(self.cov1_inpath), cov1_df.shape))

        return eqtl_df, alleles_df, cov1_df


if __name__ == '__main__':
    m = main()
    m.start()
