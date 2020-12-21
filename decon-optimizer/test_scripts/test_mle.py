#!/usr/bin/env python3

"""
File:         test_mle.py
Created:      2020/11/16
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
import argparse
import pickle
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import fmin

# Local application imports.

# Metadata
__program__ = "Test MLE"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


"""
Syntax:
./test_mle.py -p ../pickle.pkl
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.pickle_path = getattr(arguments, 'pickle')
        self.sample = getattr(arguments, 'sample')

        if not self.pickle_path.endswith(".pkl"):
            print("Data file should be a pickle file.")
            exit()

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.X = None
        self.y = None
        self.mu = None
        self.sd = None

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-p",
                            "--pickle",
                            type=str,
                            required=True,
                            help="The path to the pickle file.")
        parser.add_argument("-s",
                            "--sample",
                            type=str,
                            default="HRA_01267",
                            help="The sample to analyze.")

        return parser.parse_args()

    def start(self):
        print("### Step1 ###")
        print("Loading pickle data")
        content = self.load_pickle(self.pickle_path)
        genotype = content[0]
        cell_fractions = content[1]
        expression = content[2]

        print("### Step2 ###")
        print("Construct model.")
        self.X, self.y = self.construct_model_matrix(genotype, cell_fractions, expression)
        print("\tX:")
        print(self.X)
        print("")
        print("\ty:")
        print(self.y)
        print("")
        self.mu = self.y.mean()
        self.sd = self.y.std()

        print("### Step2 ###")
        print("Find max fast.")
        original_cf = cell_fractions.loc[self.sample]
        optimum = fmin(self.likelihood_optimization,
                       x0=np.array([original_cf]),
                       args=(self.sample, ),
                       disp=0)
        print(optimum)

        print("### Step3 ###")
        print("Calculate liklihood for different cell fractions.")
        data = []
        for cell_fraction in list(np.arange(-8, 8, 0.1)):
            data.append(self.get_log_likelihood(sample=self.sample, value=cell_fraction))

        results = pd.DataFrame(data, index=list(np.arange(-8, 8, 0.1)))
        results.reset_index(inplace=True, drop=False)
        results.columns = ["cf", "ll"]
        print(results)

        print("### Step4 ###")
        print("Plot.")
        self.plot_single(df=results,
                         x="cf",
                         y="ll",
                         xlabel="cell fraction (z-score)",
                         ylabel="log likelihood",
                         title=self.sample,
                         start_ct_frac=original_cf)

    @staticmethod
    def load_pickle(filename):
        content = {}
        i = 0
        with open(filename, "rb") as f:
            while True:
                try:
                    content[i] = pickle.load(f)
                    i += 1
                except EOFError:
                    break

        return content

    def construct_model_matrix(self, genotype, cell_fractions, y):
        # Construct the X matrix.
        X = pd.concat([genotype, cell_fractions], axis=1)
        X.columns = ["genotype", "cell_fractions"]
        X.insert(loc=0, column='intercept', value=1)
        X["genotype_x_cell_fractions"] = X["genotype"] * X["cell_fractions"]

        # Create a mask for filtering missing values.
        mask = X["genotype"] == -1

        return X.loc[~mask, :], y.loc[~mask]

    def likelihood_optimization(self, value=None, sample=None):
        return -1 * self.get_log_likelihood(sample=sample, value=value)

    def get_log_likelihood(self, sample=None, value=None):
        X = self.X
        if sample is not None and value is not None:
            X = self.change_sample_cell_fraction(sample=sample,
                                                 value=value)

        y_hat = self.create_model(X, self.y)
        ll = self.calculate_log_likelihood(y_hat)

        return ll

    def change_sample_cell_fraction(self, sample, value):
        X = self.X.copy()
        X.loc[sample, "cell_fractions"] = value
        X.loc[sample, "genotype_x_cell_fractions"] = X.loc[sample, "genotype"] * value

        return X

    def create_model(self, X, y):
        # Perform the Ordinary least squares fit.
        ols = sm.OLS(y.values, X)
        try:
            ols_result = ols.fit()
            return ols_result.predict(X)
        except np.linalg.LinAlgError as e:
            print("\tError in OLS: {}".format(e))
            return np.nan

    def calculate_log_likelihood(self, y_hat):
        return -np.sum(stats.norm.logpdf(y_hat, self.mu, self.sd))

    def plot_single(self, df, x="x", y="y", xlabel="", ylabel="", title="", start_ct_frac=None):
        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.lineplot(data=df,
                     x=x,
                     y=y,
                     ax=ax)

        best_estimate = df.loc[df[y].argmax(), x]
        ax.axvline(best_estimate, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

        start_ct_frac_string = ""
        if start_ct_frac is not None:
            ax.axvline(start_ct_frac, ls='--', color="#228B22", alpha=0.3, zorder=-1)
            start_ct_frac_string = "starting cell type fraction: {:.2f}  ".format(start_ct_frac)

        ax.text(0.5, 1.06,
                title,
                fontsize=12, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                "{}MLE estimate: {:.2f}".format(start_ct_frac_string, best_estimate),
                fontsize=10, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_xlabel(xlabel,
                      fontsize=8,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=8,
                      fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}_single_eQTL_mle_estimates.png".format(title)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Pickle path: {}".format(self.pickle_path))
        print("  > Sample: {}".format(self.sample))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
