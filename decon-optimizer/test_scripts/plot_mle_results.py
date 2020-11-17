#!/usr/bin/env python3

"""
File:         plot_mle_results.py
Created:      2020/11/17
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
import os

# Third party imports.
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Plot MLE Results"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.start_ct_frac = getattr(arguments, 'start_ct_frac')
        self.sample = getattr(arguments, 'sample')

        if not self.data_path.endswith(".pkl"):
            print("Data file should be a pickle file.")
            exit()

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data matrix (pkl).")
        parser.add_argument("-ct",
                            "--start_ct_frac",
                            type=float,
                            default=None,
                            help="The starting cell type fraction. "
                                 "Default: none.")
        parser.add_argument("-s",
                            "--sample",
                            type=str,
                            default='sample',
                            help="The starting cell type fraction. "
                                 "Default: 'sample'")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data.")
        df = pd.read_pickle(self.data_path)
        print("\tData results data frame: {}".format(df.shape))

        self.plot_multiple(df)

        cell_frac_results = df.sum(axis=0).to_frame()
        cell_frac_results.reset_index(drop=False, inplace=True)
        cell_frac_results.columns = ["x", "y"]
        self.plot_single("HRA_01267", cell_frac_results)

    def plot_multiple(self, df, nrows=4, ncols=4):
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                 figsize=(8 * nrows, 6 * ncols),
                                 sharex='all')

        row_index = 0
        col_index = 0
        for i, (index, row) in enumerate(df.iterrows()):
            ax = axes[row_index, col_index]

            subset = row.to_frame()
            subset.reset_index(drop=False, inplace=True)
            subset.columns = ["x", "y"]
            self.plot(subset, ax, index)

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

            if i == ((nrows * ncols) - 1):
                break

        fig.savefig(os.path.join(self.outdir, "{}_eqtl_level_mle_estimates.png".format(self.sample)))
        plt.close()

    def plot_single(self, title, df):
        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        self.plot(df, ax, title)

        fig.savefig(os.path.join(self.outdir, "{}_mle_estimates.png".format(title)))
        plt.close()

    def plot(self, df, ax, title):
        sns.lineplot(data=df,
                     x="x",
                     y="y",
                     ax=ax)

        best_estimate = df.loc[df['y'].argmin(), 'x']
        ax.axvline(best_estimate, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

        start_ct_frac_string = ""
        if self.start_ct_frac is not None:
            ax.axvline(self.start_ct_frac, ls='--', color="#228B22", alpha=0.3, zorder=-1)
            start_ct_frac_string = "starting cell type fraction: {}  ".format(self.start_ct_frac)

        ax.text(0.5, 1.06,
                title,
                fontsize=12, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                "{}MLE estimate: {:.2f}".format(start_ct_frac_string, best_estimate),
                fontsize=10, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_xlabel('cell fraction (z-score)',
                      fontsize=16,
                      fontweight='bold')
        ax.set_ylabel('log likelihood',
                      fontsize=16,
                      fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > Decon file: {}".format(self.data_path))
        print("  > Starting cell type fraction: {}".format(self.start_ct_frac))
        print("  > Sample name: {}".format(self.sample))
        print("  > Outpath {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
