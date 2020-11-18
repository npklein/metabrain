"""
File:         mle.py
Created:      2020/11/17
Last Changed: 2020/11/18
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
from scipy import stats
import statsmodels.api as sm

# Local application imports.


class MaximumLiklihoodEstimator:
    def __init__(self, genotype, cell_fractions, expression, log):
        # Safe arguments.
        self.X, self.y = self.construct_model_matrix(genotype,
                                                     cell_fractions,
                                                     expression)
        self.samples = list(self.X.index)
        self.log = log

    def construct_model_matrix(self, genotype, cell_fractions, y):
        # Construct the X matrix.
        X = genotype.T.merge(cell_fractions.to_frame(), left_index=True, right_index=True)
        X.columns = ["genotype", "cell_fractions"]
        X.insert(loc=0, column='intercept', value=1)
        X["genotype_x_cell_fractions"] = X["genotype"] * X["cell_fractions"]

        # Create a mask for filtering missing values.
        mask = X["genotype"] == -1

        return X.loc[~mask, :], y.loc[~mask]

    def get_n(self):
        return self.X.shape[0]

    def get_maximum_log_likelihood(self, sample=None, cell_fraction=None):
        X = self.X
        if sample is not None and cell_fraction is not None:
            if sample not in self.samples:
                return np.nan
            X = self.change_sample_cell_fraction(sample=sample,
                                                 cell_fraction=cell_fraction)

        y_hat = self.create_model(X, self.y)
        mll = self.calculate_maximum_log_likelihood(self.y, y_hat)
        # print("\t\tsample={}\tcell_fraction={}\tmll={}".format(sample, cell_fraction, mll))

        return mll

    def create_model(self, X, y):
        # Perform the Ordinary least squares fit.
        ols = sm.OLS(y.values, X)
        try:
            ols_result = ols.fit()
            return ols_result.predict(X)
        except np.linalg.LinAlgError as e:
            self.log.error("\tError in OLS: {}".format(e))
            return np.nan

    @staticmethod
    def calculate_maximum_log_likelihood(y, y_hat):
        mu = y.mean()
        sd = y.std()

        return np.log(stats.norm.pdf(y_hat, mu, sd)).sum()

    def change_sample_cell_fraction(self, sample, cell_fraction):
        X = self.X.copy()
        X.loc[sample, "cell_fractions"] = cell_fraction
        X.loc[sample, "genotype_x_cell_fractions"] = X.loc[sample, "genotype"] * cell_fraction

        return X




