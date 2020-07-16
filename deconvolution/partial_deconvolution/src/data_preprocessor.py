"""
File:         data_preprocessor.py
Created:      2020/06/29
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

# Third party imports.
import numpy as np

# Local application imports.


class DataPreprocessor:
    def __init__(self, settings, raw_signature, raw_expression):
        self.min_expr = settings.get_min_expr()
        self.normalize = settings.get_normalize()
        self.zscore = settings.get_zscore()
        self.log2 = settings.get_log2()
        self.raw_signature = raw_signature
        self.raw_expression = raw_expression

        self.shape_diff = None
        self.signature = None
        self.sign_shift = False
        self.expression = None
        self.expr_shift = False

    def work(self):
        # Filter uninformative genes from the signature matrix.
        sign_df, self.shape_diff = self.filter(self.raw_signature,
                                               self.min_expr)

        # Subset and reorder.
        sign_df, expr_df = self.subset(sign_df, self.raw_expression)

        # Transform.
        if self.normalize:
            sign_df = self.perform_normalize(sign_df, self.normalize)
        if self.log2:
            sign_df = self.perform_log2_transform(sign_df)
        if self.zscore:
            sign_df = self.perform_zscore_transform(sign_df)

        # Shift the data to be positive.
        print("Shifting data to be positive")
        if sign_df.values.min() < 0:
            sign_df = self.perform_shift(sign_df)
            self.sign_shift = True

        if expr_df.values.min() < 0:
            expr_df = self.perform_shift(expr_df)
            self.expr_shift = True

        # Save.
        self.signature = sign_df
        self.expression = expr_df

    @staticmethod
    def filter(df, cutoff):
        print("Filtering uninformative genes from the signature matrix")
        pre_shape = df.shape
        df = df.loc[(df.std(axis=1) != 0) & (df.max(axis=1) > cutoff), :]
        drop_shape = (pre_shape[0] - df.shape[0], pre_shape[1] - df.shape[1])
        return df, drop_shape

    @staticmethod
    def subset(raw_signature, raw_expression):
        print("Subsetting and reodering the input matrices")
        overlap = np.intersect1d(raw_signature.index, raw_expression.index)
        sign_df = raw_signature.loc[overlap, :]
        expr_df = raw_expression.loc[overlap, :]

        sign_df.index.name = "-"
        expr_df.index.name = "-"

        if not sign_df.index.equals(expr_df.index):
            print("Invalid order")
            exit()

        return sign_df, expr_df

    @staticmethod
    def perform_normalize(df, orient):
        print("Performing normalization:")

        if (orient == "rows") or (orient == "both"):
            print("\tper row")
            df = df.divide(df.sum(axis=1), axis=0)
        if (orient == "columns") or (orient == "both"):
            print("\tper column")
            df = df.divide(df.sum(axis=0), axis=1)
        return df

    @staticmethod
    def perform_zscore_transform(df):
        print("Performing z-score transformation")
        return df.subtract(df.mean(axis=1), axis=0).divide(df.std(axis=1), axis=0)

    @staticmethod
    def perform_log2_transform(df):
        print("Performing log2 transformation")
        df = df + abs(df.values.min()) + 1
        return np.log2(df)

    @staticmethod
    def perform_shift(df):
        return df + abs(df.values.min())

    def get_shape_diff(self):
        return self.shape_diff

    def get_signature(self):
        return self.signature

    def get_expression(self):
        return self.expression

    def get_sign_shift(self):
        return self.sign_shift

    def get_expr_shift(self):
        return self.expr_shift

    def print_info(self):
        print("Signature matrix:\t{} rows and {} columns".format(self.signature.shape[0], self.signature.shape[1]))
        print("Expression matrix:\t{} rows and {} columns".format(self.expression.shape[0], self.expression.shape[1]))
        print("Dropped:\t{} rows and {} columns".format(self.shape_diff[0], self.shape_diff[1]))
        print("Shifted:\tsignature = {}\texpression = {}".format(self.sign_shift, self.expr_shift))
        print("")
