"""
File:         storage.py
Created:      2020/10/14
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
from .container import Container


class Storage:
    def __init__(self, tech_covs, covs):
        self.tech_covs = tech_covs
        self.covs = covs

        # Initialize containers.
        self.tech_cov_container = Container(tech_covs)
        self.cov_container = Container(covs)

        # Initialize variables.
        self.error = False

    def add_row(self, eqtl_index, genotype_name):
        self.tech_cov_container.add_row(eqtl_index, genotype_name)
        self.cov_container.add_row(eqtl_index, genotype_name)

    def store_row(self):
        if not self.error:
            self.tech_cov_container.store_row()
            self.cov_container.store_row()
        else:
            print("Row not saved due to error.")

    def add_value(self, cov_name, order_id, category, type, value):
        if cov_name in self.tech_covs:
            container = self.tech_cov_container
        elif cov_name in self.covs:
            container = self.cov_container
        else:
            print("Unrecognised covariate name.")
            self.error = True
            return

        if category == "snp" and type == "coeff":
            container.add_snp_coefficient(order_id, value)
        elif category == "snp" and type == "std_err":
            container.add_snp_std_error(order_id, value)
        elif category == "inter" and type == "coeff":
            container.add_inter_coefficient(order_id, value)
        elif category == "inter" and type == "std_err":
            container.add_inter_std_error(order_id, value)
        elif category == "pvalue":
            container.add_pvalue(order_id, value)
        else:
            print("Unrecognised value characteristics: "
                  "category = {}\ttype = {}.".format(category, type))
            self.error = True
            return

    def has_error(self):
        return self.error

    def set_error(self):
        self.error = True

    def get_tech_cov_container(self):
        return self.tech_cov_container

    def get_cov_container(self):
        return self.cov_container

    def print_info(self):
        print("Storage info:")
        print("\t{} Technical covariates".format(len(self.tech_covs)))
        print("\t{} Covariates of interest".format(len(self.covs)))
        print("")
