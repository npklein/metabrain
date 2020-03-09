#!/usr/bin/env python3

"""
File:         plot_z_scores.py
Created:      2020/02/27
Last Changed: 2020/03/03
Author:       M.Vochteloo

Copyright (C) 2019 M.Vochteloo

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
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.


# Metadata.
__program__ = "Plot Z-scores"
__author__ = "M. Vochteloo"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class Main:
    """
    Main class of the program.
    """

    def __init__(self, regr_file):
        """
        Initializer method for the main class.

        :param eqtl_file: string, the file containing the regression data.
        """
        self.regr_file = regr_file
        self.outdir = os.path.join(os.getcwd(), 'output', 'plots')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        sns.set()

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Load the genotype matrix file.
        print("Loading regression matrix.")
        data = pd.read_csv(self.regr_file, sep="\t", header=0, nrows=nrows)
        print("\tShape: {}".format(data.shape))

        # Calculate the correlation.
        print("Calculating correlation.")
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            data["overal_z_score"],
            data["z_score_estimate"])

        # # Subset the weird cases.
        # tmp = data.sort_values(by=['z_score_estimate'], ascending=False)
        # print(tmp.iloc[0:10, :])
        # exit()

        # Add colors.
        data["hue"] = "#2C7BB6"
        data.loc[data["flipped"], "hue"] = "#000000"

        # Print.
        print("Info:")
        print("\tSlope: {:.4f}\tIntercept: {:.4f}\t"
              "Correlation coefficient: {:.4f}\t"
              "P-value: {:.2e}\tStandard error: {:.2e}\t".format(slope,
                                                                 intercept,
                                                                 r_value,
                                                                 p_value,
                                                                 std_err
                                                                 ))

        # Plot.
        print("Plotting z-scores.")
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        g = sns.regplot(x="z_score_estimate", y="overal_z_score",
                        data=data,
                        scatter_kws={'facecolors': data['hue'],
                                     'edgecolor': data['hue'],
                                     'alpha': 0.5},
                        line_kws={"color": "#D7191C"},
                        )
        g.set_title('Z-score regression (Pearson r = {:.2f}, '
                    'p = {:.2e})'.format(r_value, p_value))
        g.set_ylabel('Overal Z-score',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('Estimated Z-score (slope / std_err)',
                     fontsize=8,
                     fontweight='bold')
        g.axhline(0, ls='--', color="#000000", alpha=0.3, zorder=-1)
        g.axvline(0, ls='--', color="#000000", alpha=0.3, zorder=-1)
        #plt.show()
        fig.savefig(os.path.join(self.outdir, "z_score_regression.png"))

        # Plot the p-value distribution.
        print("Plotting p-value distribution.")
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        g = sns.distplot(data[["p_value"]])
        g.set_title('P-value Distribution')
        g.set_ylabel('Frequency')
        g.set_xlabel('P-value')
        #plt.show()
        fig.savefig(os.path.join(self.outdir, "pvalue_distribution.png"))

        # Plot the beta distribution.
        print("Plotting slope distribution.")
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        g = sns.distplot(data[["slope"]])
        g.set_title('Slope Distribution')
        g.set_ylabel('Frequency')
        g.set_xlabel('Slope (-1 > x < 1)')
        plt.xlim(-1, 1)
        #plt.show()
        fig.savefig(os.path.join(self.outdir, "slope_distribution.png"))


if __name__ == "__main__":
    # Define main variables.
    REGRESSION = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output", "2019-11-06-FreezeTwoDotOne",
                              "2020-03-03-zscore-comparison",
                              "step1-regression-matrix", "output",
                              "regression_table.txt.gz")

    # Start the program.
    MAIN = Main(regr_file=REGRESSION)
    MAIN.start()
