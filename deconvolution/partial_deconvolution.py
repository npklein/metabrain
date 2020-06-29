#!/usr/bin/env python3

"""
File:         partial_deconvolution.py
Created:      2020/06/29
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

# Third party imports.

# Local application imports.
from partial_deconvolution.src.main import Main
from partial_deconvolution.src.settings import Settings
from partial_deconvolution.src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Partial Deconvolution"
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
if __name__ == '__main__':
    # Get the command line arguments.
    CLA = CommandLineArguments(program=__program__,
                               version=__version__,
                               description=__description__)
    SETTINGS = Settings(data_path=CLA.get_argument("data"),
                        signature_path=CLA.get_argument("signature"),
                        translate_path=CLA.get_argument("translate"),
                        sample_path=CLA.get_argument("sample"),
                        cohort=CLA.get_argument("cohort"),
                        min_expr=CLA.get_argument("min_expr"),
                        normalize=CLA.get_argument("normalize"),
                        zscore=CLA.get_argument("zscore"),
                        log2=CLA.get_argument("log2"),
                        decon_method=CLA.get_argument("decon_method"),
                        sum_to_one=CLA.get_argument("sum_to_one"),
                        extension=CLA.get_argument("extension")
                        )

    # Start the program.
    PROGRAM = Main(settings=SETTINGS, outdir=CLA.get_argument("outdir"))
    PROGRAM.start()
