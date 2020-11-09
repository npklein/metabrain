"""
File:         cmd_line_arguments.py
Created:      2020/11/09
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
import argparse

# Third party imports.

# Local application imports.


class CommandLineArguments:
    def __init__(self, program, version, description):
        # Safe variables.
        self.program = program
        self.version = version
        self.description = description

        # Get the arguments.
        parser = self.create_argument_parser()
        self.arguments = parser.parse_args()
        self.print_arguments()

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=self.program,
                                         description=self.description)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(self.program,
                                                   self.version),
                            help="show program's version number and exit")
        parser.add_argument("-i",
                            "--input_dir_path",
                            type=str,
                            default=None,
                            help="The path to the cell type specific expression"
                                 "matrices input directory. "
                                 "Default: 'output'.")
        parser.add_argument("-s",
                            "--input_suffix",
                            type=str,
                            default="tsv",
                            help="The suffix of the cell type specific expression"
                                 "matrices. Default: 'tsv'.")
        parser.add_argument("-r",
                            "--reference_profile_path",
                            type=str,
                            required=True,
                            help="The cell type specific reference profile "
                                 "for deconvolution.")
        parser.add_argument("-g",
                            "--gene_info_path",
                            type=str,
                            required=True,
                            help="The path to the gene info file.")
        parser.add_argument("-n",
                            "--n_samples",
                            type=int,
                            default=1000,
                            help="The number of artifical bulk samples to "
                                 "create. Default: 1000.")
        parser.add_argument("-c",
                            "--combine",
                            nargs="*",
                            type=str,
                            default=None,
                            help="Option to combine cell types. Syntax = "
                                 "cellType1-cellType2-Name.")

        return parser

    def print_arguments(self):
        for arg in vars(self.arguments):
            print("Input argument '{}' "
                  "has value '{}'.".format(arg, getattr(self.arguments, arg)))

    def get_argument(self, arg_key):
        if self.arguments is not None and arg_key in self.arguments:
            value = getattr(self.arguments, arg_key)
        else:
            value = None

        return value

    def get_all_arguments(self):
        return self.arguments
