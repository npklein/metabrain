#!/usr/bin/env python3

"""
File:         plot_batch_correction.py
Created:      2020/12/21
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
__program__ = "Plot batch Correction"
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

        if not self.data_path.endswith(".pkl"):
            print("Data file should be a pickle file.")
            exit()

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.info_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-03-09.brain.phenotypes.txt"
        self.palette = {
            "AMP-AD": "#A9AAA9",
            "Braineac": "#E7A023",
            "Brainseq": "#CC79A9",
            "CMC": "#EEE643",
            "GTEx": "#1A9F74",
            "NABEC": "#767833",
            "TargetALS": "#DDCC78",
            "ENA": "#D56228",
            "PsychEncode": "#5AB4E5"
        }

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

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data.")
        df = pd.read_pickle(self.data_path).to_frame()
        print("\tData results data frame: {}".format(df.shape))
        df.columns = ["cell fraction"]
        print(df)

        info_df = pd.read_csv(self.info_path, sep="\t", header=0,
                              index_col=4, low_memory=False)
        print("\tInformation data frame: {} "
              "with shape: {}".format(os.path.basename(self.info_path),
                                      info_df.shape))
        print(info_df)

        print("### Step2 ###")
        print("Add cohort.")
        sample_metacohort_dict = dict(zip(info_df.index, info_df["MetaCohort"]))
        df["MetaCohort"] = df.index.map(sample_metacohort_dict)
        df["Color"] = df["MetaCohort"].map(self.palette)
        print(df)

        print("### Step3 ###")
        print("Plotting.")
        self.plot(df, self.palette, self.outdir)

    @staticmethod
    def plot(df, palette, outdir):
        # Plot.
        fig, ax = plt.subplots(figsize=(15, 8))
        sns.despine(fig=fig, ax=ax)

        sns.set()
        sns.set_style("ticks")

        sns.scatterplot(x=df.index,
                        y=df["cell fraction"],
                        hue=df["MetaCohort"],
                        palette=palette,
                        linewidth=0,
                        legend="brief",
                        ax=ax)

        ax.axhline(df["cell fraction"].mean(), ls='--', color="#b22222", zorder=-1)

        ax.set_ylabel("cell fraction",
                      fontsize=10,
                      fontweight='bold')
        ax.set_xlabel("",
                      fontsize=10,
                      fontweight='bold')

        plt.setp(ax.get_legend().get_texts(), fontsize='2')
        plt.setp(ax.get_legend().get_title(), fontsize='4')
        ax.legend(bbox_to_anchor=(0.5, -0.45), loc='lower center',
                  ncol=5)

        fig.savefig(os.path.join(outdir, "batch_correction.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Pickle file: {}".format(self.data_path))
        print("  > Outpath {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
