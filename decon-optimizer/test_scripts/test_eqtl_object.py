#!/usr/bin/env python3

"""
File:         test_eqtl_object.py
Created:      2020/11/21
Last Changed: 2020/11/23
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
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Local application imports.

# Metadata
__program__ = "Test eQTL Object"
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
./test_eqtl_object.py -p ../pickle.pkl
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

        print("### Step3 ###")
        print("Get the coef representation.")
        # start_time = int(time.time())
        # for i in range(327):
        #     coefs = self.get_coef_representation()
        # rt_min, rt_sec = divmod(int(time.time()) - start_time, 60)
        # print("\t\tfinished in {} minute(s) and "
        #       "{} second(s)".format(int(rt_min),
        #                             int(rt_sec)))

        coefs = self.get_coef_representation()
        vertex = self.calculate_vertex(coefs)
        if vertex is not None:
            print("\tVertex: ({:.4f}, {:.4f})".format(vertex[0], vertex[1]))
        pol_func = np.poly1d(coefs)

        print("### Step4 ###")
        print("Find max fast.")
        original_cf = cell_fractions.loc[self.sample]
        optimum = self.optimize_cell_fraction(self.sample)
        print(optimum)

        print("### Step5 ###")
        print("Calculate liklihood for different cell fractions.")
        data = []
        x_values = list(np.arange(-8, 8, 0.1))
        for cell_fraction in x_values:
            real = self.get_log_likelihood(sample=self.sample, value=cell_fraction)
            pred = pol_func(cell_fraction)
            data.append([real, pred])

        results = pd.DataFrame(data, index=x_values, columns=["ll", "func"])
        results.index.name = "cf"
        results.reset_index(inplace=True, drop=False)
        results = results.melt(id_vars=["cf"])
        print(results)

        print("### Step6 ###")
        print("Plot.")
        self.line_plot(df=results,
                       x="cf",
                       y="value",
                       hue="variable",
                       xlabel="cell fraction (z-score)",
                       ylabel="log likelihood",
                       title=self.sample,
                       start_ct_frac=original_cf)

        print("### Step7 ###")
        print("Plot the interaction model.")
        inter_df = self.X.copy()
        inter_df["hue"] = inter_df["genotype"].round(0)
        inter_df = inter_df[["cell_fractions", "hue"]].merge(self.y, left_index=True, right_index=True)
        print(inter_df)

        sample_df = inter_df.loc[[self.sample], ["cell_fractions", self.y.name]]
        sample_df.columns = ["cf", "expr"]
        sample_df.index = ["original"]
        sample_df.loc["calculation"] = [float(vertex[0]), sample_df.at["original", "expr"]]
        sample_df.loc["optimization"] = [float(optimum), sample_df.at["original", "expr"]]
        sample_df.reset_index(drop=False, inplace=True)
        print(sample_df)

        self.inter_plot(df=inter_df,
                        sample_df=sample_df,
                        x="cell_fractions",
                        y=self.y.name,
                        hue="hue",
                        xlabel="cell fraction (z-score)",
                        ylabel=self.y.name,
                        title="{} interaction eQTL".format(self.sample),
                        subtitle="genotype: {}".format(self.X.at[self.sample, "genotype"]))

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

    @staticmethod
    def construct_model_matrix(genotype, cell_fractions, y):
        # Construct the X matrix.
        X = pd.concat([genotype, cell_fractions], axis=1)
        X.columns = ["genotype", "cell_fractions"]
        X.insert(loc=0, column='intercept', value=1)
        X["genotype_x_cell_fractions"] = X["genotype"] * X["cell_fractions"]

        # Create a mask for filtering missing values.
        mask = X["genotype"] == -1

        return X.loc[~mask, :], y.loc[~mask]

    def optimize_cell_fraction(self, sample):
        res = minimize(self.likelihood_optimization,
                       x0=np.array([self.X.at[sample, "cell_fractions"]]),
                       method="Powell",
                       args=(sample,),
                       tol=0.001,
                       options={"disp": False})
        return res.x

    def get_derivate_func(self):
        return np.polyder(np.poly1d(self.get_coef_representation()))

    def likelihood_optimization(self, value=None, sample=None):
        return -1 * self.get_log_likelihood(sample=sample, value=value)

    def get_coef_representation(self):
        # evaluate the parabula at degree+1 points. Since it is quadratic,
        # pick three points.
        x_values = [-8, self.X.at[self.sample, "cell_fractions"], 8]
        y_values = []
        for value in x_values:
            y_values.append(self.get_log_likelihood(sample=self.sample,
                                                    value=value))

        return np.polyfit(x_values, y_values, 2)

    @staticmethod
    def calculate_vertex(coefs):
        if len(coefs) != 3:
            return None
        return (-coefs[1] / (2 * coefs[0])), (((4 * coefs[0] * coefs[2]) - (coefs[1] * coefs[1])) / (4 * coefs[0]))

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


    def line_plot(self, df, x="x", y="y", hue=None, xlabel="", ylabel="", title="", start_ct_frac=None):
        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.lineplot(data=df,
                     x=x,
                     y=y,
                     hue=hue,
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

    def inter_plot(self, df, sample_df, x="x", y="y", hue=None, xlabel="", ylabel="", title="", subtitle=""):
        if hue is None:
            return
        if len(set(df[hue].unique()).symmetric_difference({0, 1, 2})) > 0:
            return

        colors = {0.0: "#E69F00",
                  1.0: "#0072B2",
                  2.0: "#D55E00"}

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        label_pos = {0.0: 0.94, 1.0: 0.90, 2.0: 0.86}
        for i, group in enumerate(df[hue].unique()):
            subset = df.loc[df[hue] == group, :].copy()
            color = colors[group]

            coef_str = "NA"
            p_str = "NA"
            if len(subset.index) > 1:
                # Regression.
                coef, p = stats.spearmanr(subset[x], subset[y])
                coef_str = "{:.2f}".format(coef)
                p_str = "p = {:.2e}".format(p)

                # Plot.
                sns.regplot(x=x, y=y, data=subset,
                            scatter_kws={'facecolors': color,
                                         'linewidth': 0,
                                         'alpha': 0.15},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            # Add the text.
            ax.annotate(
                '{}: r = {} [{}]'.format(group, coef_str, p_str),
                # xy=(0.03, 0.94 - ((i / 100) * 4)),
                xy=(0.03, label_pos[group]),
                xycoords=ax.transAxes,
                color=color,
                alpha=0.75,
                fontsize=12,
                fontweight='bold')

        # Plot the sample start point.
        sns.scatterplot(x="cf", y="expr", hue="index", data=sample_df)

        ax.text(0.5, 1.06,
                title,
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                subtitle,
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}.png".format(title.replace(" ", "_"))))
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
