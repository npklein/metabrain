"""
File:         utilities.py
Created:      2020/03/12
Last Changed:
Author(s):    M.Vochteloo

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
from pathlib import Path
import glob
import os

# Third party imports.

# Local application imports.


def get_project_root_dir():
    """
    Method to get the project root dir.

    :return: string, project root directory.
    """
    path = Path(__file__).parent.parent
    return str(path)


def check_file_exists(file_path):
    """
    Method to check if the file path exists.

    :param file_path: string, the complete file path.
    :return: bool, True if file exists; False if not.
    """
    return os.path.exists(file_path)


def get_extension(file_path):
    """
    Method to get the extension of a file_path.

    :param file_path: string, the complete file path.
    :return: string, the extension of the file.
    """
    return os.path.splitext(file_path)[1][1:]


def get_file_name(file_path):
    """
    Method to get the base name (filename inc. extension) of a file_path

    :param file_path: string, the complete file path.
    :return: string, the filename.
    """
    return os.path.basename(file_path)


def prepare_output_dir(out_dir):
    """
    Method to create a directory if it does not exist.

    :param out_dir: string, a directory.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def p_value_to_symbol(p_value):
    """
    Method for translating a p-value to a symbol.

    :param p_value: float, the input p-value.
    :return: string, the symbol corresponding to the p-value.
    """
    output = ""
    try:
        if p_value > 0.05:
            output = "NS"
        elif 0.05 >= p_value > 0.01:
            output = "*"
        elif 0.01 >= p_value > 0.001:
            output = "**"
        elif 0.001 >= p_value > 0.0001:
            output = "***"
        elif p_value <= 0.0001:
            output = "****"
    except TypeError:
        pass

    return output
