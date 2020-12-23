"""
File:         eqtl_object.py
Created:      2020/11/17
Last Changed: 2020/12/23
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
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

# Local application imports.


class EQTLObject:
    def __init__(self, genotype, cell_fractions, expression, log):
        # Safe arguments.
        self.X, self.y = self.construct_model_matrix(genotype,
                                                     cell_fractions,
                                                     expression)
        self.log = log

        # Calculate distribution properties.
        self.mu = self.y.mean()
        self.sd = self.y.std()

        # Set standard values.
        self.x0, self.x2 = -8, 8

    @staticmethod
    def construct_model_matrix(genotype, cell_fractions, y):
        X = pd.concat([genotype, cell_fractions], axis=1)
        X.columns = ["genotype", "cell_fractions"]
        X.insert(loc=0, column='intercept', value=1)
        X["genotype_x_cell_fractions"] = X["genotype"] * X["cell_fractions"]

        # Create a mask for filtering missing values.
        mask = X["genotype"] == -1

        return X.loc[~mask, :], y.loc[~mask]

    def contains_sample(self, sample):
        return sample in self.X.index

    def calculate_optimal_cell_fraction(self, sample):
        coefs = self.get_ll_coef_representation(sample)
        return self.calculate_vertex_x_coordinate(a=coefs[0], b=coefs[1])

    def get_ll_coef_representation(self, sample):
        # evaluate the parabula at degree+1 points. Since it is quadratic,
        # pick three points.
        x_values = [self.x0, self.X.at[sample, "cell_fractions"], self.x2]
        y_values = []
        for value in x_values:
            y_values.append(self.get_log_likelihood(sample=sample,
                                                    value=value))

        return np.polyfit(x_values, y_values, 2)

    @staticmethod
    def calculate_vertex_x_coordinate(a, b):
        return -b / (2 * a)

    def optimize_cell_fraction(self, sample):
        res = minimize(self.likelihood_optimization,
                       x0=np.array([self.X.at[sample, "cell_fractions"]]),
                       method="Powell",
                       args=(sample,),
                       tol=0.001,
                       options={"disp": False})
        return res.x

    def likelihood_optimization(self, value=None, sample=None):
        return -1 * self.get_log_likelihood(sample=sample, value=value)

    def get_log_likelihood(self, sample=None, value=None):
        X = self.X
        if sample is not None and value is not None and value != X.at[sample, "cell_fractions"]:
            X = self.change_sample_cell_fraction(sample=sample,
                                                 value=value)

        y_hat = self.fit_and_predict_mlr_model(X, self.y)
        ll = self.calculate_log_likelihood(y_hat)

        return ll

    def change_sample_cell_fraction(self, sample, value):
        X = self.X.copy()
        X.at[sample, "cell_fractions"] = value
        X.at[sample, "genotype_x_cell_fractions"] = X.at[sample, "genotype"] * value

        return X

    @staticmethod
    def fit_and_predict_mlr_model(X, y):
        return X.dot(np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y))

    def calculate_log_likelihood(self, y_hat):
        return -np.sum(stats.norm.logpdf(y_hat, self.mu, self.sd))





