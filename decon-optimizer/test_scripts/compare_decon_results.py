#!/usr/bin/env python3

"""
File:         compare_decon_results.py
Created:      2020/12/28
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
from statsmodels.stats import multitest
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Compare Decon Results"
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
./compare_decon_results.py -od /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -ad /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-10-decon-optimizer/2020-12-28-decon-eQTL/cis/cortex/decon_out/deconvolutionResults.csv -ocf /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt -acf /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-10-decon-optimizer/data/cell_fractions.txt -ocn CellMapNNLS_Neuron -acn CellMapNNLS_Neuron_ocf 
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.ori_decon_path = getattr(arguments, 'ori_decon')
        self.alt_decon_path = getattr(arguments, 'alt_decon')
        self.ori_cf_path = getattr(arguments, 'ori_cell_fractions')
        self.alt_cf_path = getattr(arguments, 'alt_cell_fractions')
        self.ori_colname = getattr(arguments, 'ori_colname')
        self.alt_colname = getattr(arguments, 'alt_colname')

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
        parser.add_argument("-od",
                            "--ori_decon",
                            type=str,
                            required=True,
                            help="The path to the original deconvolution "
                                 "matrix.")
        parser.add_argument("-ad",
                            "--alt_decon",
                            type=str,
                            required=True,
                            help="The path to the alternative deconvolution "
                                 "matrix.")
        parser.add_argument("-ocf",
                            "--ori_cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the original cell fractions "
                                 "matrix.")
        parser.add_argument("-acf",
                            "--alt_cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the alternative cell fractions "
                                 "matrix.")
        parser.add_argument("-ocn",
                            "--ori_colname",
                            type=str,
                            required=True,
                            help="The name of the original results column.")
        parser.add_argument("-acn",
                            "--alt_colname",
                            type=str,
                            required=True,
                            help="The name of the alterative results column.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading original data.")
        ori_decon_df, ori_cf_df = self.load_data(self.ori_decon_path, self.ori_cf_path, self.ori_colname)
        print("Loading alternative data.")
        alt_decon_df, alt_cf_df = self.load_data(self.alt_decon_path, self.alt_cf_path, self.alt_colname)

        if self.ori_colname == self.alt_colname:
            ori_decon_df.columns = ["ori_{}".format(x) for x in ori_decon_df.columns]
            ori_cf_df.columns = ["ori_{}".format(x) for x in ori_cf_df.columns]
            ori_decon_df.columns = ["alt_{}".format(x) for x in ori_decon_df.columns]
            alt_cf_df.columns = ["alt_{}".format(x) for x in alt_cf_df.columns]

            self.ori_colname = "ori_{}".format(self.ori_colname)
            self.alt_colname = "alt_{}".format(self.alt_colname)

        print("### Step2 ###")
        print("Calculate BH-FDR.")
        ori_decon_df["{}_FDR".format(self.ori_colname)] = multitest.multipletests(ori_decon_df.loc[:, "{}_pvalue".format(self.ori_colname)], method='fdr_bh')[1]
        alt_decon_df["{}_FDR".format(self.alt_colname)] = multitest.multipletests(alt_decon_df.loc[:, "{}_pvalue".format(self.alt_colname)], method='fdr_bh')[1]

        print("### Step3 ###")
        print("Merge.")
        decon_df = ori_decon_df.merge(alt_decon_df, left_index=True, right_index=True)
        cf_df = ori_cf_df.merge(alt_cf_df, left_index=True, right_index=True)
        print(decon_df)
        print(cf_df)

        print("### Step4 ###")
        print("Calculate significant hits.")
        counts = decon_df[decon_df < 0.05].count()
        for index in counts.index:
            print("\t'{}' has {} values < 0.05.".format(index, counts[index]))

        print("### Step5 ###")
        print("Plotting cell fractions.")
        self.scatterplot(data=cf_df,
                         x=self.ori_colname,
                         y=self.alt_colname,
                         xlabel="original cell fraction",
                         ylabel="optimized cell fraction",
                         title="cell fraction comparison",
                         filename="cell_fraction_comparison")

        print("### Step6 ###")
        print("Plotting FDR values.")
        self.scatterplot(data=decon_df,
                         x="{}_FDR".format(self.ori_colname),
                         y="{}_FDR".format(self.alt_colname),
                         xlabel="original cf - FDR",
                         ylabel="optimized cf - FDR",
                         title="deconvolution FDR comparison",
                         filename="decon_fdr_comparison")

        print("### Step7 ###")
        print("Plotting significant hits.")
        sign_decon_df = decon_df.loc[decon_df["{}_FDR".format(self.ori_colname)] < 0.01, :].copy()
        self.scatterplot(data=sign_decon_df,
                         x="{}_FDR".format(self.ori_colname),
                         y="{}_FDR".format(self.alt_colname),
                         xlabel="original cf - FDR",
                         ylabel="optimized cf - FDR",
                         title="deconvolution FDR comparison",
                         filename="sign_decon_fdr_comparison")

    @staticmethod
    def load_data(decon_path, cf_path, colname):
        decon_df = pd.read_csv(decon_path,
                               sep="\t",
                               header=0,
                               index_col=0)
        print("\tDeconvolution data frame: {}".format(decon_df.shape))

        cf_df = pd.read_csv(cf_path,
                            sep="\t",
                            header=0,
                            index_col=0)
        print("\tCell fractions data frame: {}".format(cf_df.shape))

        return decon_df[["{}_pvalue".format(colname)]], cf_df[[colname]]

    def scatterplot(self, data, x="x", y="y", xlabel="", ylabel="", title="",
                    filename="image"):
        # Plot.
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")

        # ax.axvline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)
        # ax.axhline(0, ls='--', color="#000000", alpha=0.15, zorder=-1)

        g = sns.scatterplot(data=data,
                            x=x,
                            y=y,
                            legend=False,
                            alpha=0.5)

        g.set_title(title)
        g.set_ylabel(ylabel,
                     fontsize=10,
                     fontweight='bold')
        g.set_xlabel(xlabel,
                     fontsize=10,
                     fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Original decon path: {}".format(self.ori_decon_path))
        print("  > Alternative decon path: {}".format(self.alt_decon_path))
        print("  > Original cell type fraction path: {}".format(self.ori_cf_path))
        print("  > Alternative cell type fraction path: {}".format(self.alt_cf_path))
        print("  > Original column name: {}".format(self.ori_colname))
        print("  > Alternative column name: {}".format(self.alt_colname))
        print("  > Outpath {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
