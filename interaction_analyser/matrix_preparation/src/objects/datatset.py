"""
File:         dataset.py
Created:      2020/03/12
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

# Third party imports.

# Local application imports.
from src.df_utilities import load_dataframe


class Dataset:
    def __init__(self, geno_path, alleles_path, expr_path, cov_path, nrows):
        self.geno_path = geno_path
        self.alleles_path = alleles_path
        self.expr_path = expr_path
        self.cov_path = cov_path
        self.nrows = nrows

        self.geno_df = None
        self.alleles_df = None
        self.expr_df = None
        self.cov_df = None

    def get_genotype(self):
        if self.geno_df is None:
            self.geno_df = load_dataframe(self.geno_path,
                                          header=0,
                                          index_col=0,
                                          nrows=self.nrows)
        return self.geno_df

    def get_alleles(self):
        if self.alleles_df is None:
            self.alleles_df = load_dataframe(self.alleles_path,
                                             header=0,
                                             index_col=0,
                                             nrows=self.nrows)
        return self.geno_df

    def get_expression(self):
        if self.expr_df is None:
            self.expr_df = load_dataframe(self.expr_path,
                                          header=0,
                                          index_col=0,
                                          nrows=self.nrows)
        return self.geno_df

    def get_covariates(self):
        if self.cov_df is None:
            self.cov_df = load_dataframe(self.cov_path,
                                         header=0,
                                         index_col=0,
                                         nrows=self.nrows)
        return self.geno_df
