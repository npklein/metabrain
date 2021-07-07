#!/usr/bin/env python3

"""
File:         cellmap_reference_profile.py
Created:      2021/07/06
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
import gzip
import os

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Cellmap Reference Profile"
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
./cellmap_reference_profile.py -d ../data/CellMap_brain_CNS7_log2CPM.txt -o CellMap_brain_CNS7_avgCPM.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.outfile = getattr(arguments, 'outfile')

        self.cell_type_dict = {
            "Oligodendrocytes": "Oligodendrocyte",
            "Endothelial": "EndothelialCell",
            "Pericyte": "Pericytes",
            "Astrocytes": "Astrocyte"
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
                            help="The path to the Cellmap reference profile")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the output file")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        df = self.load_file(self.data_path, header=0, index_col=0)

        print("Preprocessing header")
        info = []
        for col in df.columns:
            info.append(col.split("|"))
        info_df = pd.DataFrame(info, columns=["cell type", "dataset", "num", "iter"])
        print(info_df)

        print("Step 1: inversing log2 transform.")
        df = df.pow(2)

        print("Step 2: take average per cell type")
        avg_df = None
        for cell_type in info_df["cell type"].unique():
            mask = (info_df["cell type"] == cell_type).to_numpy()
            ct_avg_expr = df.loc[:, mask].mean(axis=1).to_frame()
            adj_cell_type = cell_type
            if adj_cell_type in self.cell_type_dict:
                adj_cell_type = self.cell_type_dict[cell_type]
            ct_avg_expr.columns = ["CellMap_{}".format(adj_cell_type)]

            if avg_df is None:
                avg_df = ct_avg_expr
            else:
                avg_df = pd.concat([avg_df, ct_avg_expr], axis=1)

        print("Step 3: prepare output matrix")
        avg_df.insert(0, "GeneSymbol", avg_df.index)

        print("Saving file")
        self.save_file(avg_df,
                       outpath=os.path.join(os.path.dirname(self.data_path), self.outfile),
                       index=False)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Data path: {}".format(self.data_path))
        print("  > Output filename: {}".format(self.outfile))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
