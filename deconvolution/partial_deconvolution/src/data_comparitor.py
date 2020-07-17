"""
File:         data_comparitor.py
Created:      2020/06/30
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
import pandas as pd

# Local application imports.


class DataComparitor:
    def __init__(self, settings, deconvolution, ground_truth):
        self.deconvolution = deconvolution
        self.ground_truth = ground_truth

        self.comparison = None
        self.n_samples = None
        self.rss = None

    def work(self):
        if self.ground_truth is None:
            return

        decon_df = self.deconvolution.copy()
        truth_df = self.ground_truth.copy()

        overlap = np.intersect1d(decon_df.index, truth_df.index)
        decon_df = decon_df.loc[overlap, :]
        truth_df = truth_df.loc[overlap, :]

        decon_df.index.name = "-"
        truth_df.index.name = "-"

        if not decon_df.index.equals(truth_df.index):
            print("Invalid order")
            exit()

        cell_types = np.intersect1d(decon_df.columns, truth_df.columns)
        data = []
        for i, sample in enumerate(overlap):
            decon = decon_df.iloc[i, :]
            truth = truth_df.iloc[i, :]
            for cell in cell_types:
                data.append([sample, decon[cell], truth[cell], cell])
        df = pd.DataFrame(data, columns=['-', 'x', 'y', 'hue'])
        df.set_index('-', inplace=True)

        self.comparison = df
        self.n_samples = len(overlap)
        self.rss = sum((df['x'] - df['y'])**2)

    def get_comparison(self):
        return self.comparison

    def get_n_samples(self):
        return self.n_samples

    def get_rss(self):
        return self.rss

    def print_info(self):
        if self.comparison is not None:
            print("Comparison matrix:\t{} rows and {} columns".format(self.comparison.shape[0], self.comparison.shape[1]))
        if self.n_samples is not None:
            print("N samples:\t{}".format(self.n_samples))
        if self.rss is not None:
            print("Residuals sum of squares:\t{:.2f}".format(self.rss))
        print("")
