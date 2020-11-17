"""
File:         cell_type_object.py
Created:      2020/11/16
Last Changed: 2020/11/17
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
import time

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats

# Local application imports.
from .mle import MaximumLiklihoodEstimator


class CellType:
    def __init__(self, cell_type, eqtl_df, geno_df, alleles_df, expr_df,
                 ct_fractions, cohort_sample_dict, log):
        # Safe arguments.
        self.cell_type = cell_type
        self.eqtl_df = eqtl_df
        self.eqtl_df.reset_index(drop=True, inplace=True)
        self.n = self.eqtl_df.shape[0]
        self.geno_df = geno_df
        self.alleles_df = alleles_df
        self.cohort_sample_dict = cohort_sample_dict
        self.log = log

        # Force normal on expression and cell count matrices.
        self.expr_df = self.force_normal_per_cohort(expr_df)
        self.ct_frac_df = self.force_normal(ct_fractions, value_axis=0)

        # Validate.
        self.validate()

        # Declare empty variables.
        self.mle_objects = self.create_mle_objects()

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

        if not self.ct_frac_df.index.equals(self.geno_df.columns):
            self.log.error("The genotype file columns do not match the cell "
                           "type fractions file.")
            exit()

        if not self.ct_frac_df.index.equals(self.expr_df.columns):
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

    def create_mle_objects(self):
        mle_objects = {}
        for index, row in self.eqtl_df.iterrows():
            key = row["SNPName"] + "_" + row["ProbeName"]
            mle_objects[key] = MaximumLiklihoodEstimator(genotype=self.geno_df.iloc[[index], :].copy(),
                                                         cell_fractions=self.ct_frac_df,
                                                         expression=self.expr_df.iloc[index, :].copy(),
                                                         log=self.log)

        return mle_objects

    def test_all_cell_fractions(self, sample):
        self.log.info("\tTesting all possible cell type fractions for sample "
                      "'{}'".format(sample))
        start_time = int(time.time())

        all_data = []
        indices = []
        columns = list(np.arange(-4, 4, 0.1))
        for i, (eqtl_name, mle_object) in enumerate(self.mle_objects.items()):
            if (i % 25 == 0) or (i == (self.n - 1)):
                self.log.info("\t\tProcessing eQTL {} / {} [{:.2f}%]".format(i, self.n - 1, (100/(self.n - 1)) * i))
            row_data = []
            for cell_fraction in columns:
                mll = mle_object.get_maximum_log_likelihood(sample=sample,
                                                            cell_fraction=cell_fraction)
                row_data.append(mll)

            all_data.append(row_data)
            indices.append(eqtl_name)

        results = pd.DataFrame(all_data, index=indices, columns=columns)
        results.dropna(axis=0, inplace=True)
        cell_frac_results = results.sum(axis=0)
        print(cell_frac_results)

        self.log.info("\tReal fraction: {:.4f}\tmll fraction: {:.4f}".format(self.ct_frac_df.loc[sample][0],
                                                                             cell_frac_results.idxmax()))

        run_time_min, run_time_sec = divmod(int(time.time()) - start_time, 60)
        self.log.info("\tFinished in {} minute(s) and "
                      "{} second(s)".format(int(run_time_min),
                                            int(run_time_sec)))

        return results



