"""
File:         dataset.py
Created:      2020/03/16
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
import itertools

# Third party imports.

# Local application imports.
from src.df_utilities import load_dataframe


class Dataset:
    def __init__(self, settings):
        self.eqtl_inpath = settings.get_setting("eqtl_datafile")
        self.geno_inpath = settings.get_setting("genotype_datafile")
        self.alleles_inpath = settings.get_setting("alleles_datafile")
        self.expr_inpath = settings.get_setting("expression_datafile")
        self.cov_inpath = settings.get_setting("covariate_datafile")
        self.inter_inpath = settings.get_setting("interaction_datafile")
        self.nrows = settings.get_setting("nrows")

        # Declare empty variables.
        self.eqtl_df = None
        self.geno_df = None
        self.alleles_df = None
        self.expr_df = None
        self.cov_df = None
        self.inter_df = None

    def get_eqtl_df(self):
        if self.eqtl_df is None:
            self.eqtl_df = load_dataframe(inpath=self.eqtl_inpath,
                                          header=0,
                                          index_col=False,
                                          nrows=self.nrows)
            self.eqtl_df.index = self.eqtl_df["SNPName"]
            self.eqtl_df.index.name = "-"
            self.validate()
        return self.eqtl_df

    def get_geno_df(self):
        if self.geno_df is None:
            self.geno_df = load_dataframe(inpath=self.geno_inpath,
                                          header=0,
                                          index_col=0,
                                          nrows=self.nrows)
            self.validate()
        return self.geno_df

    def get_alleles_df(self):
        if self.alleles_df is None:
            self.alleles_df = load_dataframe(inpath=self.alleles_inpath,
                                             header=0,
                                             index_col=0,
                                             nrows=self.nrows)
            self.validate()
        return self.alleles_df

    def get_expr_df(self):
        if self.expr_df is None:
            self.expr_df = load_dataframe(inpath=self.expr_inpath,
                                          header=0,
                                          index_col=0,
                                          nrows=self.nrows)
            self.validate()
        return self.expr_df

    def get_cov_df(self):
        if self.cov_df is None:
            self.cov_df = load_dataframe(inpath=self.cov_inpath,
                                         header=0,
                                         index_col=0)
            self.validate()
        return self.cov_df

    def get_inter_df(self):
        if self.inter_df is None:
            self.inter_df = load_dataframe(inpath=self.inter_inpath,
                                           header=0,
                                           index_col=0)
            self.validate()
        return self.inter_df

    def validate(self):
        dfs = [self.eqtl_df, self.geno_df, self.alleles_df, self.expr_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.index.identical(b.index):
                print("Order of eQTLs is not identical (1).")
                exit()

        dfs = [self.geno_df, self.expr_df, self.cov_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.columns.identical(b.columns):
                print("Order of samples are not identical.")
                exit()

        if self.cov_df is not None and self.inter_df is not None and \
                not self.cov_df.index.identical(self.inter_df.index):
            print("Order of covariates is not identical.")
            exit()

        if self.inter_df is not None:
            subset = self.inter_df.iloc[:, :self.nrows].copy()
            for i, colname in enumerate(subset.columns):
                for df in [self.eqtl_df, self.geno_df, self.alleles_df,
                           self.expr_df]:
                    if not colname.startswith(df.index[i]):
                        print("Order of eQTLs is not identical (2).")
                        exit()
            del subset

        print("\tValid.")
