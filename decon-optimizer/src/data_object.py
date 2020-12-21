"""
File:         data_object.py
Created:      2020/11/16
Last Changed: 2020/12/21
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
import os

# Third party imports.
import pandas as pd
from statsmodels.stats import multitest

# Local application imports.


class Data:
    def __init__(self, eqtl_path, genotype_path, alleles_path, expression_path,
                 cell_fractions_path, decon_path,  sample_annotation_path,
                 sample_id, cohort_id, log):
        # Safe arguments.
        self.eqtl_path = eqtl_path
        self.geno_path = genotype_path
        self.alle_path = alleles_path
        self.expr_path = expression_path
        self.frac_path = cell_fractions_path
        self.deco_path = decon_path
        self.saan_path = sample_annotation_path
        self.sample_id = sample_id
        self.cohort_id = cohort_id
        self.log = log

        # Set empty variables.
        self.eqtl_df = None
        self.geno_df = None
        self.alle_df = None
        self.expr_df = None
        self.frac_df = None
        self.deco_df = None
        self.deco_fdr_df = None
        self.saan_df = None

    def get_eqtl_df(self, skiprows=None, nrows=None):
        if self.eqtl_df is None:
            self.eqtl_df = self.load_dataframe(self.eqtl_path,
                                               header=0,
                                               index_col=None,
                                               skiprows=skiprows,
                                               nrows=nrows)

        return self.eqtl_df

    def get_geno_df(self, skiprows=None, nrows=None):
        if self.geno_df is None:
            self.geno_df = self.load_dataframe(self.geno_path,
                                               header=0,
                                               index_col=0,
                                               skiprows=skiprows,
                                               nrows=nrows)

        return self.geno_df

    def get_alle_df(self, skiprows=None, nrows=None):
        if self.alle_df is None:
            self.alle_df = self.load_dataframe(self.alle_path,
                                               header=0,
                                               index_col=0,
                                               skiprows=skiprows,
                                               nrows=nrows)

        return self.alle_df

    def get_expr_df(self, skiprows=None, nrows=None):
        if self.expr_df is None:
            self.expr_df = self.load_dataframe(self.expr_path,
                                               header=0,
                                               index_col=0,
                                               skiprows=skiprows,
                                               nrows=nrows)

        return self.expr_df

    def get_frac_df(self, skiprows=None, nrows=None):
        if self.frac_df is None:
            self.frac_df = self.load_dataframe(self.frac_path,
                                               header=0,
                                               index_col=0,
                                               skiprows=skiprows,
                                               nrows=nrows)

        return self.frac_df

    def get_deco_df(self, skiprows=None, nrows=None):
        if self.deco_df is None:
            deco_df = self.load_dataframe(self.deco_path,
                                          header=0,
                                          index_col=0,
                                          skiprows=skiprows,
                                          nrows=nrows)

            # Split index into ProbeName and SNPName
            probe_names = []
            snp_names = []
            for index in deco_df.index:
                probe_names.append(index.split("_")[0])
                snp_names.append("_".join(index.split("_")[1:]))
            deco_df["ProbeName"] = probe_names
            deco_df["SNPName"] = snp_names
            deco_df.reset_index(drop=True, inplace=True)

            self.deco_df = deco_df

        return self.deco_df

    def get_deco_fdr_df(self):
        if self.deco_fdr_df is None:
            if self.deco_df is None:
                _ = self.get_deco_df()

            data = []
            indices = []
            for colname, values in self.deco_df.T.iterrows():
                if colname in ["ProbeName", "SNPName"]:
                    data.append(values)
                    indices.append(colname)
                if colname.endswith("_pvalue"):
                    data.append(multitest.multipletests(values, method='fdr_bh')[1])
                    indices.append(colname.replace("_pvalue", ""))
            self.deco_fdr_df = pd.DataFrame(data, index=indices, columns=self.deco_df.index).T

        return self.deco_fdr_df

    def load_dataframe(self, inpath, header, index_col, sep="\t", low_memory=True,
                       nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        self.log.info("\tLoaded dataframe: {} "
                      "with shape: {}".format(os.path.basename(inpath),
                                              df.shape))
        return df

    def get_cohort_to_sample_dict(self):
        sample_cohort_dict = self.get_sample_annotation_dict(key=self.sample_id,
                                                             value=self.cohort_id)

        return self.reverse_merge_dict(sample_cohort_dict)

    @staticmethod
    def reverse_merge_dict(dict):
        out_dict = {}
        seen_keys = set()
        for key, value in dict.items():
            if key in seen_keys:
                print("Key {} has muiltiple values.".format(key))
            seen_keys.add(key)

            if value in out_dict.keys():
                keys = out_dict[value]
                keys.append(key)
                out_dict[value] = keys
            else:
                out_dict[value] = [key]

        return out_dict

    def get_sample_annotation_dict(self, key, value):
        if self.saan_df is None:
            self.saan_df = self.load_dataframe(self.saan_path,
                                               header=0,
                                               index_col=None,
                                               low_memory=False)

        return dict(zip(self.saan_df.loc[:, key], self.saan_df.loc[:, value]))

    def print_arguments(self):
        self.log.info("Data Arguments:")
        self.log.info("  > eQTL input path: {}".format(self.eqtl_path))
        self.log.info("  > Genotype input path: {}".format(self.geno_path))
        self.log.info("  > Alleles input path: {}".format(self.alle_path))
        self.log.info("  > Expression input path: {}".format(self.expr_path))
        self.log.info("  > Cell fractions input path: {}".format(self.frac_path))
        self.log.info("  > Decon input path: {}".format(self.deco_path))
        self.log.info("  > Sample annotation path: {}".format(self.saan_path))
        self.log.info("  > Sample ID: {}".format(self.sample_id))
        self.log.info("  > Cohort ID: {}".format(self.cohort_id))
        self.log.info("")
