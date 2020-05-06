"""
File:         container.py
Created:      2020/05/06
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


class Container:
    def __init__(self, colnames):
        # Initialize the result list.
        self.pvalues_buffer = [[-1, "-"] + colnames]
        self.tvalues_buffer = [[-1, "-"] + colnames]
        self.perm_pvalues_buffer = []
        self.perm_tvalues_buffer = []

        # Initialize a new row.
        self.pvalues = None
        self.tvalues = None
        self.perm_pvalues = None
        self.perm_tvalues = None

    def add_row(self, eqtl_index, genotype_name):
        self.pvalues = [eqtl_index, genotype_name]
        self.tvalues = [eqtl_index, genotype_name]

        self.perm_pvalues = []
        self.perm_tvalues = []

    def store_row(self):
        self.pvalues_buffer.append(self.pvalues)
        self.tvalues_buffer.append(self.tvalues)
        self.perm_pvalues_buffer.extend(self.perm_pvalues)
        self.perm_tvalues_buffer.extend(self.perm_tvalues)

        self.pvalues = None
        self.tvalues = None
        self.perm_pvalues = None
        self.perm_tvalues = None

    def add_tvalue(self, order_id, value):
        if order_id == 0:
            self.tvalues.append(value)
        else:
            self.perm_tvalues.append(value)

    def add_pvalue(self, order_id, value):
        if order_id == 0:
            self.pvalues.append(value)
        else:
            self.perm_pvalues.append(value)

    def get_pvalues(self):
        return self.pvalues_buffer

    def get_tvalues(self):
        return self.tvalues_buffer

    def get_perm_pvalues(self):
        return self.perm_pvalues_buffer

    def get_perm_tvalues(self):
        return self.perm_tvalues_buffer
