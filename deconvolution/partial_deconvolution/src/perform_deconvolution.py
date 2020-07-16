"""
File:         perform_deconvolution.py
Created:      2020/06/29
Last Changed: 2020/06/30
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
import os

# Third party imports.
from scipy.optimize import nnls
import pandas as pd

# Local application imports.


class PerformDeconvolution:
    def __init__(self, settings, signature, expression):
        self.decon_method = settings.get_decon_method()
        self.sum_to_one = settings.get_sum_to_one()
        self.outdir = settings.get_output_path()
        self.signature = signature
        self.expression = expression

        self.deconvolution = None
        self.residuals = None

    def work(self):
        print("Performing deconvolution using '{}'".format(self.decon_method))
        decon_function = None
        if self.decon_method == "NNLS":
            decon_function = self.nnls
        else:
            print("Unexpected deconvolution method.")
            exit()

        decon_data = []
        residuals_data = []
        for index, sample in self.expression.T.iterrows():
            proportions, rnorm = decon_function(self.signature, sample)
            decon_data.append(proportions)
            residuals_data.append(rnorm)

        deconvolution = pd.DataFrame(decon_data,
                                     index=self.expression.columns,
                                     columns=self.signature.columns)
        residuals = pd.Series(residuals_data,
                              index=self.expression.columns)

        # Make weights sum up to one.
        if self.sum_to_one:
            deconvolution = self.perform_sum_to_one(deconvolution)

        # Save.
        self.deconvolution = deconvolution
        self.residuals = residuals
        self.save()

    @staticmethod
    def nnls(A, b):
        return nnls(A, b)

    @staticmethod
    def perform_sum_to_one(X):
        print("Making weights sum-to-one")
        return X.divide(X.sum(axis=1), axis=0)

    def get_deconvolution(self):
        return self.deconvolution

    def get_residuals(self):
        return self.residuals

    def get_avg_residuals(self):
        return self.residuals.mean()

    def get_info_per_celltype(self):
        means = self.deconvolution.mean(axis=0).to_dict()
        stds = self.deconvolution.std(axis=0).to_dict()

        info = {}
        for key in means.keys():
            info[key] = (means[key], stds[key])
        return info

    def save(self):
        self.deconvolution.to_csv(os.path.join(self.outdir, 'deconvolution.txt.gz'),
                                  compression="gzip",
                                  sep="\t",
                                  header=True,
                                  index=True)
        print("\tsaved dataframe: deconvolution.txt.gz "
              "with shape: {}".format(self.deconvolution.shape))

    def print_info(self):
        print("Average residuals: {:.2f}".format(self.get_avg_residuals()))
        print("Average weights per celltype:")
        for celltype, (mean, std) in self.get_info_per_celltype().items():
            print("\t{:20s}: mean = {:.2f} std = {:.2f}".format(celltype, mean, std))
        print("")