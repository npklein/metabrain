"""
File:         cell_type_object.py
Created:      2020/11/16
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
import pandas as pd
from scipy import stats

# Local application imports.


class CellType:
    def __init__(self, cell_type, eqtl_df, geno_df, alleles_df, expr_df,
                 ct_fractions, cohort_sample_dict, log):
        # Safe arguments.
        self.cell_type = cell_type
        self.eqtl_df = eqtl_df
        self.geno_df = geno_df
        self.alleles_df = alleles_df
        self.cohort_sample_dict = cohort_sample_dict
        self.log = log

        # Force normal on expression and cell count matrices.
        self.expr_df = self.force_normal_per_cohort(expr_df)
        self.ct_fractions = self.force_normal(ct_fractions, value_axis=0)

        # Validate.
        self.validate()

    def validate(self):
        snp_reference = self.eqtl_df["SNPName"].copy()
        snp_reference.rename("-", inplace=True)

        probe_reference = self.eqtl_df["ProbeName"].copy()
        probe_reference.rename("-", inplace=True)

        if not pd.Series(self.geno_df.index,
                         index=snp_reference.index,
                         name="-").equals(snp_reference):
            self.log.error("The genotype file indices do not match the eQTL "
                           "file.")
            exit()

        if not pd.Series(self.alleles_df.index,
                         index=snp_reference.index,
                         name="-").equals(snp_reference):
            self.log.error("The genotype alleles file indices do not match the "
                           "eQTL file.")
            exit()

        if not pd.Series(self.expr_df.index,
                         index=snp_reference.index,
                         name="-").equals(probe_reference):
            self.log.error("The expression file indices do not match the eQTL "
                           "file.")
            exit()

        if not self.ct_fractions.index.equals(self.geno_df.columns):
            self.log.error("The genotype file columns do not match the cell "
                           "type fractions file.")
            exit()

        if not self.ct_fractions.index.equals(self.expr_df.columns):
            self.log.error("The expressiom file columns do not match the cell "
                           "type fractions file.")
            exit()

    def force_normal_per_cohort(self, df, sample_axis=0):
        work_df = df.copy()

        value_axis = 1
        if sample_axis == 0:
            pass
        elif sample_axis == 1:
            value_axis = 0
            work_df = work_df.T
        else:
            self.log.error("Unexpected axis in force normal per cohort "
                           "function.")
            exit()

        sample_order = df.columns

        combined_df = None
        for cohort, cohort_samples in self.cohort_sample_dict.items():
            mask = []
            for sample in sample_order:
                if sample not in cohort_samples:
                    mask.append(False)
                else:
                    mask.append(True)

            cohort_df = work_df.loc[:, mask]
            if cohort_df.shape[1] == 0:
                continue

            normal_df = self.force_normal(cohort_df, value_axis)

            if combined_df is None:
                combined_df = normal_df
            else:
                combined_df = pd.concat([combined_df, normal_df], axis=1)

        if sample_axis == 1:
            combined_df = combined_df.T

        return combined_df.loc[:, sample_order]

    def force_normal(self, df, value_axis=1):
        work_df = df.copy()

        if value_axis == 0:
            work_df = work_df.T
        elif value_axis == 1:
            pass
        else:
            self.log.error("Unexpected axis in force normal function.")
            exit()

        data = []
        for index, row in work_df.iterrows():
            data.append(self.force_normal_series(row))

        normal_df = pd.DataFrame(data, index=work_df.index, columns=work_df.columns)

        if value_axis == 0:
            normal_df = normal_df.T

        return normal_df

    @staticmethod
    def force_normal_series(s):
        return stats.norm.ppf((s.rank(ascending=True) - 0.5) / s.size)

