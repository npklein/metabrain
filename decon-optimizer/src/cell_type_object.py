"""
File:         cell_type_object.py
Created:      2020/11/16
Last Changed: 2020/12/23
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
from .eqtl_object import EQTLObject
from .utilities import force_normal_series
from .multiprocessor import MultiProcessor
from scipy.optimize import minimize


class CellType:
    def __init__(self, cell_type, eqtl_df, alleles_df, indices, cohorts,
                 sample_order, cores, max_end_time, log):
        # Safe arguments.
        self.cell_type = cell_type
        self.indices = list(eqtl_df.index)
        self.eqtl_df = eqtl_df
        self.n = eqtl_df.shape[0]
        self.alleles_df = alleles_df
        self.indices = indices
        self.cohorts = cohorts
        self.sample_order = sample_order
        self.cores = cores
        self.max_end_time = max_end_time
        self.log = log

        # Create other arguments.
        self.geno_df = self.get_combined_matrix(df="geno")
        self.expr_df = self.get_combined_matrix(df="expr", type="normal")
        #self.cf_s = self.get_combined_matrix(df="cf_df", type="normal")
        self.cf_s = force_normal_series(self.get_combined_matrix(df="cf_df"), as_series=True)

        # Create MLE objects.
        self.eqtl_objects = self.create_eqtl_objects()

        # Create empty variables.
        self.mp = None

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

    def create_eqtl_objects(self):
        eqtl_objects = {}
        for index, row in self.eqtl_df.iterrows():
            key = row["SNPName"] + "_" + row["ProbeName"]
            eqtl_objects[key] = EQTLObject(genotype=self.geno_df.iloc[index, :].copy(),
                                           cell_fractions=self.cf_s,
                                           expression=self.expr_df.iloc[index, :].copy(),
                                           log=self.log)

        return eqtl_objects

    def test_cell_fraction_range(self, poi=None):
        self.log.info("\tTesting all possible cell type fractions")
        results = self.perform_analysis(work_func=self.work_function_test_cell_fraction_range,
                                        poi=poi)

        return self.combine_dict_of_series_into_df(results)

    def calculate_optimal_cell_fraction(self, poi=None):
        self.log.info("\tCalculating optimal log likelihood cell fraction")
        results = self.perform_analysis(work_func=self.work_function_calculate_optimal_cell_fraction,
                                        poi=poi)
        return pd.Series(results)

    def optimize_cell_fraction_per_eqtl(self, poi=None):
        self.log.info("\tOptimizing maximum log likelihood cell fraction per eQTL")
        results = self.perform_analysis(work_func=self.work_function_optimize_cell_fraction_per_eqtl,
                                        poi=poi)

        return self.combine_dict_of_series_into_df(results)

    def optimize_cell_fraction(self, poi=None):
        self.log.info("\tOptimizing maximum log likelihood cell fractions")
        results = self.perform_analysis(work_func=self.work_function_optimize_cell_fraction,
                                        poi=poi)
        return pd.Series(results)

    @staticmethod
    def combine_dict_of_series_into_df(dict):
        columns = []
        series = []
        for sample, s in dict.items():
            columns.append(sample)
            series.append(s)

        df = pd.concat(series, axis=1)
        df.columns = columns
        df.index.name = "-"
        return df

    def perform_analysis(self, work_func, poi):
        start_time = int(time.time())

        if isinstance(poi, str):
            (sample, output) = work_func(sample=poi)
            results = {sample: output}
        elif isinstance(poi, list) or isinstance(poi, set) or isinstance(poi, pd.Series):
            if self.mp is None:
                self.mp = MultiProcessor(samples=self.sample_order,
                                         cores=self.cores,
                                         max_end_time=self.max_end_time,
                                         log=self.log)
            results = self.mp.process(work_func=work_func, samples=poi)
        else:
            if self.mp is None:
                self.mp = MultiProcessor(samples=self.sample_order,
                                         cores=self.cores,
                                         max_end_time=self.max_end_time,
                                         log=self.log)
            results = self.mp.process(work_func=work_func)

        rt_min, rt_sec = divmod(int(time.time()) - start_time, 60)
        self.log.info("\t\tfinished in {} minute(s) and "
                      "{} second(s)".format(int(rt_min),
                                            int(rt_sec)))

        return results

    def work_function_test_cell_fraction_range(self, sample):
        sample_sll = []
        sample_indices = []
        for cf in list(np.arange(-8, 8, 0.1)):
            sample_sll.append(self.get_summed_log_likelihood(sample=sample,
                                                             value=cf))
            sample_indices.append(cf)

        return pd.Series(sample_sll, index=sample_indices)

    def work_function_calculate_optimal_cell_fraction(self, sample):
        a = 0
        b = 0
        for eqtl_object in self.eqtl_objects.values():
            if eqtl_object.contains_sample(sample):
                coefs = eqtl_object.get_ll_coef_representation(sample)
                a += coefs[0]
                b += coefs[1]

        return self.calculate_vertex_x_coordinate(a, b)

    @staticmethod
    def calculate_vertex_x_coordinate(a, b):
        return -b / (2 * a)

    def work_function_optimize_cell_fraction(self, sample):
        res = minimize(self.likelihood_optimization,
                       x0=np.array([self.cf_s.at[sample]]),
                       method="Powell",
                       args=(sample,),
                       tol=0.001,
                       options={"disp": False})
        return res.x

    def work_function_optimize_cell_fraction_per_eqtl(self, sample):
        ll_values = []
        indices = []
        for eqtl_name, eqtl_object in self.eqtl_objects.items():
            if eqtl_object.contains_sample(sample):
                ll_values.append(eqtl_object.optimize_cell_fraction(sample))
                indices.append(eqtl_name)
        return pd.Series(ll_values, index=indices)

    def likelihood_optimization(self, value, sample):
        return -1 * self.get_summed_log_likelihood(value=value, sample=sample)

    def get_summed_log_likelihood(self, value, sample):
        total_ll = 0
        for eqtl_object in self.eqtl_objects.values():
            if eqtl_object.contains_sample(sample):
                total_ll += eqtl_object.get_log_likelihood(sample=sample,
                                                           value=value)
        return total_ll

    def print_info(self):
        self.log.info("\tCelltype: {}".format(self.cell_type))
        self.log.info("\t > N-eQTLs: {}".format(self.eqtl_df.shape[0]))
        self.log.info("\t > N-samples: {}".format(len(self.sample_order)))

