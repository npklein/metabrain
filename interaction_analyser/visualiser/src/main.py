"""
File:         main.py
Created:      2020/03/13
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
import os

# Third party imports.

# Local application imports.
from src.utilities import get_project_root_dir, prepare_output_dir
from src.local_settings import LocalSettings


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, outdir):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param force_steps: list, the names of the steps to force to redo.
        :param outdir: string, the name of the base output directory.
        """
        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(settings_file)

        # Prepare an output directory.
        self.outdir = os.path.join(get_project_root_dir(), outdir)
        prepare_output_dir(self.outdir)

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting visualiser.")
        self.print_arguments()

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("")
