#!/usr/bin/env python3

"""
File:         search_dataframe.py
Created:      2020/07/16
Last Changed: 2020/07/17
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
import argparse

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Search DataFrame"
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
        arguments = self.create_argument_parser()
        self.snps = getattr(arguments, 'snps')
        self.genes = getattr(arguments, 'genes')
        self.data_path = getattr(arguments, 'data')
        self.eqtl_path = getattr(arguments, 'eqtls')
        self.filter = getattr(arguments, 'filter')
        self.header = getattr(arguments, 'header')
        self.index_col = getattr(arguments, 'index_col')
        self.sep = getattr(arguments, 'sep')
        self.alpha = getattr(arguments, 'alpha')

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
        parser.add_argument("-s",
                            "--snps",
                            nargs="+",
                            type=str,
                            default=[],
                            help="The SNPs of interest. Default: [].")
        parser.add_argument("-g",
                            "--genes",
                            nargs="+",
                            type=str,
                            default=[],
                            help="The genes of interest. Default: [].")
        parser.add_argument("-f",
                            "--filter",
                            nargs="+",
                            type=str,
                            default=None,
                            help="The indices to filter on. Default: None.")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            help="The data to look in.")
        parser.add_argument("-e",
                            "--eqtls",
                            type=str,
                            help="The eqtls corresponding to the data.")
        parser.add_argument("--header",
                            type=int,
                            default=0,
                            help="The data header row. Default: 0.")
        parser.add_argument("--index_col",
                            type=int,
                            default=0,
                            help="The data index column. Default: 0.")
        parser.add_argument("--sep",
                            type=str,
                            default="\t",
                            help="The data column separator. Default: '\\t'.")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            default=0.05,
                            help="The significance cut-off. Default: 0.05.")

        return parser.parse_args()

    def start(self):
        print("Load the data.")
        data_df = pd.read_csv(self.data_path,
                              header=self.header,
                              index_col=self.index_col,
                              sep=self.sep)

        eqtl_df = pd.read_csv(self.eqtl_path,
                              header=self.header,
                              index_col=self.index_col,
                              sep=self.sep)

        print("Searching rows.")
        for i, (index, row) in enumerate(data_df.iterrows()):
            if self.is_hit(index):
                self.print_info(i, index, row, eqtl_df.iloc[i, :])

        print("Searching columns.")
        for i, (index, row) in enumerate(data_df.T.iterrows()):
            if self.is_hit(index):
                self.print_info(i, index, row, eqtl_df.iloc[i, :])

    def is_hit(self, index):
        for snp in self.snps:
            if snp in index:
                return True
        for gene in self.genes:
            if gene in index:
                return True
        return False

    def print_info(self, i, name, row, eqtl):
        probe_name = ""
        hgnc_name = ""
        if name.startswith(eqtl["SNPName"]):
            probe_name = eqtl["ProbeName"]
            hgnc_name = eqtl["HGNCName"]
        buffer = ["Index:\t{} - {} [{}]".format(name, probe_name, hgnc_name)]
        for index, item in row.iteritems():
            if self.filter is None:
                if item < self.alpha:
                    buffer.append("\t[{}]:{}\t{:.4f}".format(i, index, item))
            else:
                for value in self.filter:
                    if value in index:
                        if item < self.alpha:
                            buffer.append("\t[{}]:{}\t{:.4f}".format(i, index, item))
                            break
        buffer.append("")

        if len(buffer) > 2:
            for line in buffer:
                print(line)


if __name__ == '__main__':
    m = main()
    m.start()
