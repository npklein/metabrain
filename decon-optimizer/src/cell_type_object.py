"""
File:         cell_type_object.py
Created:      2020/11/16
Last Changed: 2020/11/18
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

# Local application imports.
from .mle import MaximumLiklihoodEstimator


class CellType:
    def __init__(self, cell_type, eqtl_df, alleles_df, indices, cohorts,
                 sample_order, log):
        # Safe arguments.
        self.cell_type = cell_type
        self.indices = list(eqtl_df.index)
        self.eqtl_df = eqtl_df
        self.n = eqtl_df.shape[0]
        self.alleles_df = alleles_df
        self.indices = indices
        self.cohorts = cohorts
        self.sample_order = sample_order
        self.log = log

        # Create other arguments.
        self.geno_df = self.get_combined_matrix(df="geno")
        print(self.geno_df)
        self.expr_df = self.get_combined_matrix(df="expr", type="normal")
        print(self.expr_df)
        self.cf_s = self.get_combined_matrix(df="cf_df", type="normal")
        print(self.cf_s)

        # Create MLE objects.
        self.mle_objects = self.create_mle_objects()

    def get_combined_matrix(self, df, type=None):
        combined_df = None
        for cohort_object in self.cohorts.values():
            normal_df = None
            if df == "geno" and type is None:
                normal_df = cohort_object.get_geno_df()
            elif df == "geno" and type == "normal":
                normal_df = cohort_object.get_normal_expr_df()
            elif df == "expr" and type is None:
                normal_df = cohort_object.get_expr_df()
            elif df == "expr" and type == "normal":
                normal_df = cohort_object.get_normal_expr_df()
            elif df == "cf_df" and type is None:
                normal_df = cohort_object.get_cf_df()
            elif df == "cf_df" and type == "normal":
                normal_df = cohort_object.get_normal_cf_df()
            else:
                pass

            if normal_df is None:
                self.log.error("Unexpeted input for get_combined_matrix function.")
                exit()

            if combined_df is None:
                combined_df = normal_df
            else:
                if df == "cf_df":
                    combined_df = pd.concat([combined_df, normal_df], axis=0)
                else:
                    combined_df = pd.concat([combined_df, normal_df], axis=1)

        if df == "cf_df":
            # Reorder samples.
            combined_df = combined_df.loc[self.sample_order]
            # Subset cell type.
            combined_df = combined_df.loc[:, self.cell_type]
        else:
            # Reorder samples.
            combined_df = combined_df.loc[:, self.sample_order]
            # Subset eQTLs.
            combined_df = combined_df.iloc[self.indices, :]

        return combined_df

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

        if not self.cf_s.index.equals(self.geno_df.columns):
            self.log.error("The genotype file columns do not match the cell "
                           "type fractions file.")
            exit()

        if not self.cf_s.index.equals(self.expr_df.columns):
            self.log.error("The expressiom file columns do not match the cell "
                           "type fractions file.")
            exit()

    def create_mle_objects(self):
        mle_objects = {}
        for index, row in self.eqtl_df.iterrows():
            key = row["SNPName"] + "_" + row["ProbeName"]
            mle_objects[key] = MaximumLiklihoodEstimator(genotype=self.geno_df.iloc[[index], :].copy(),
                                                         cell_fractions=self.cf_s,
                                                         expression=self.expr_df.iloc[index, :].copy(),
                                                         log=self.log)

        return mle_objects

    def test_all_cell_fractions(self, sample):
        self.log.info("\tTesting all possible cell type fractions for sample "
                      "'{}'".format(sample))
        start_time = int(time.time())

        all_data = []
        indices = []
        columns = list(np.arange(0, 1, 0.025))
        for i, (eqtl_name, mle_object) in enumerate(self.mle_objects.items()):
            if (i % 25 == 0) or (i == (self.n - 1)):
                self.log.info("\t\tProcessing eQTL {} / {} [{:.2f}%]".format(i, self.n - 1, (100/(self.n - 1)) * i))
            row_data = []
            for cell_fraction in columns:
                original_value, normalized_cell_fraction = self.calculate_new_normalized_ct_fraction(sample, cell_fraction)
                mll = mle_object.get_maximum_log_likelihood(sample=sample,
                                                            cell_fraction=normalized_cell_fraction)
                # print("Original value: {:.4f}\tcell fraction: {:.4f}\t"
                #       "normalized cell fraction: {:.4f}\t"
                #       "mll: {}".format(original_value,
                #                        cell_fraction,
                #                        normalized_cell_fraction,
                #                        mll))
                row_data.append(mll)

            all_data.append(row_data)
            indices.append(eqtl_name)

        results = pd.DataFrame(all_data, index=indices, columns=columns)
        results.dropna(axis=0, inplace=True)
        print(results)
        cell_frac_results = results.sum(axis=0)
        print(cell_frac_results)

        self.log.info("\tReal fraction: {:.4f}\tmll fraction: {:.4f}".format(self.cf_s.loc[sample],
                                                                             cell_frac_results.idxmax()))

        run_time_min, run_time_sec = divmod(int(time.time()) - start_time, 60)
        self.log.info("\tFinished in {} minute(s) and "
                      "{} second(s)".format(int(run_time_min),
                                            int(run_time_sec)))

        return results

    def calculate_new_normalized_ct_fraction(self, sample, new_value):
        original_value = None
        new_normalized_value = None
        for cohort, cohort_object in self.cohorts.items():
            if cohort_object.contains_sample(sample):
                original_value, new_normalized_value = cohort_object.calculate_new_normalized_ct_fraction(sample, self.cell_type, new_value)
                break

        return original_value, new_normalized_value

    def print_info(self):
        self.log.info("Celltype: {}".format(self.cell_type))
        self.log.info(" > N-eQTLs: {}".format(self.eqtl_df.shape[0]))
        self.log.info(" > N-samples: {}".format(self.eqtl_df.shape[1]))
        self.log.info("")

