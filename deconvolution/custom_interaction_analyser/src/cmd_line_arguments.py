"""
File:         cmd_line_arguments.py
Created:      2020/03/23
Last Changed: 2020/04/02
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
    """
    CommandLineArguments: class that processing command line arguments.
    """

    def __init__(self, program, version, description):
        """
        Initializer method.

        :param program: string, name of the program.
        :param version: float, the version of the program.
        :param description: string, description of the program.
        """
        # Safe variables.
        self.program = program
        self.version = version
        self.description = description

        # Get the arguments.
        parser = self.create_argument_parser()
        self.arguments = parser.parse_args()
        self.print()

    def create_argument_parser(self):
        """
        Function to create the command line flags.

        :return: parser object,
        """
        parser = argparse.ArgumentParser(prog=self.program,
                                         description=self.description)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(self.program,
                                                   self.version),
                            help="show program's version number and exit")
        parser.add_argument("-s",
                            "--settings",
                            type=str,
                            default="default_settings",
                            help="The settings input file (without '.json'), "
                                 "default: 'default_settings'.")
        parser.add_argument("-sr",
                            "--skip_rows",
                            type=int,
                            default=0,
                            help="The number of rows to skip in the input "
                                 "files, default: 0.")
        parser.add_argument("-ne",
                            "--n_eqtls",
                            type=int,
                            default=None,
                            help="The number of eQTLs in the input files, "
                                 "default: None (determine manually).")
        parser.add_argument("-ns",
                            "--n_samples",
                            type=int,
                            default=None,
                            help="The number of samples in the input files, "
                                 "default: None (determine manually).")
        parser.add_argument("-c",
                            "--cores",
                            type=int,
                            default=1,
                            help="The number of cores to use, "
                                 "default: 1.")
        parser.add_argument("-verbose",
                            action='store_true',
                            help="Include steps and command prints, "
                                 "default: False.")
        parser.add_argument("-combine",
                            action='store_true',
                            help="Combine the created files, alternative "
                                 "functionality. Default: False.")

        return parser

    def print(self):
        """
        Method that prints the input settings.

        :return:
        """
        for arg in vars(self.arguments):
            print("Input argument '{}' "
                  "has value '{}'.".format(arg, getattr(self.arguments, arg)))

    def get_argument(self, arg_key):
        """
        Method to return the value of a command line key.

        :param arg_key: str, key in parse command line
        :return: int/str/dict, value of arg_key. If the key did not
                 exist, return None.
        """
        if self.arguments is not None and arg_key in self.arguments:
            value = getattr(self.arguments, arg_key)
        else:
            value = None

        return value

    def get_all_arguments(self):
        """
        Method to return all command line arguments.

        :return self.arguments: dict, dictionary containing all command
                line arguments.
        """
        return self.arguments
