#!/usr/bin/env python3

"""
File:         visualise_deconvolution_matrix.py
Created:      2020/09/04
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
__program__ = "Visualise Deconvolution Matrix"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
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
        self.decon_path = getattr(arguments, 'decon')
        self.info_path = getattr(arguments, 'info')
        self.sample_id = getattr(arguments, 'sample')
        self.x_id = getattr(arguments, 'x')
        self.row_id = getattr(arguments, 'row')

        # Set variables.
        self.outdir = str(Path(__file__).parent.parent)
        self.colormap = {
            "Neuron": "#b38d84",
            "Oligodendrocyte": "#5d9166",
            "EndothelialCell": "#f2a7a7",
            "Microglia": "#e8c06f",
            "Macrophage": "#e8c06f",
            "Astrocyte": "#9b7bb8"
        }

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-d",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the deconvolution matrix")
        parser.add_argument("-i",
                            "--info",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the sample information matrix")
        parser.add_argument("-sample",
                            type=str,
                            required=True,
                            help="The sample column.")
        parser.add_argument("-x",
                            type=str,
                            required=True,
                            help="The information column for the x-axis.")
        parser.add_argument("-row",
                            type=str,
                            required=True,
                            help="The information column for the rows.")

        return parser.parse_args()

    def start(self):
        print("Loading deconvolution matrix.")
        decon_df = pd.read_csv(self.decon_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon_path),
                                      decon_df.shape))

        print("Loading information matrix.")
        info_df = pd.read_csv(self.info_path, sep="\t", header=0,
                              low_memory=False)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.info_path),
                                      info_df.shape))

        print("Preprocessing.")
        x_dict = dict(zip(info_df.loc[:, self.sample_id], info_df.loc[:, self.x_id]))
        row_dict = dict(zip(info_df.loc[:, self.sample_id], info_df.loc[:, self.row_id]))

        df = decon_df.reset_index().melt(id_vars=["index"])
        df[self.x_id] = df["index"].map(x_dict).astype(str)
        df[self.row_id] = df["index"].map(row_dict).astype(str)

        print("Visualizing.")
        sns.set(style="ticks")
        order = list(df[self.x_id].unique())
        order.sort()
        g = sns.catplot(x=self.x_id, y="value", col="variable",
                        row=self.row_id, data=df, kind="box",
                        order=order)
        [plt.setp(ax.texts, text="") for ax in g.axes.flat]
        #g.set_titles(row_template='{row_name}', col_template='{col_name}')

        for axes in g.axes.flat:
            axes.set_xticklabels(axes.get_xticklabels(),
                                 rotation=65,
                                 horizontalalignment='right')
            old_title = axes.get_title().replace(" ", "")

            row_column, row_variable = old_title.split("|")[0].split("=")
            col_column, col_variable = old_title.split("|")[1].split("=")

            subset = df[(df[row_column] == row_variable) & (df[col_column] == col_variable)]
            counts = subset[self.x_id].value_counts()

            counts_string = []
            total = 0
            for key in order:
                value = 0
                if key in counts:
                    value = counts[key]
                counts_string.append(str(value))
                total += int(value)

            axes.set_title("{} | {} | N = {}\n[{}]".format(row_variable, col_variable, str(total), ", ".join(counts_string)))

        g.savefig(os.path.join(self.outdir, "x{}_row{}_catplot.png".format(self.x_id.replace(" ", "_"), self.row_id.replace(" ", "_"))))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
