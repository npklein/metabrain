"""
File:         dataset.py
Created:      2020/03/16
Last Changed: 2020/05/12
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
import itertools
import os

# Third party imports.
import scipy.stats as stats

# Local application imports.
from general.df_utilities import load_dataframe


class Dataset:
    def __init__(self, settings, nrows):
        self.input_dir = settings.get_setting("input_dir")
        filenames = settings.get_setting("filenames")
        self.eqtl_filename = filenames["eqtl"]
        self.geno_filename = filenames["genotype"]
        self.alleles_filename = filenames["alleles"]
        self.expr_filename = filenames["expression"]
        self.cov_filename = filenames["covariates"]
        self.markers_filename = filenames["markers"]

        self.inter_input_dir = settings.get_setting("interaction_input_dir")
        inter_subdirs = settings.get_setting("interaction_input_subfolders")
        self.inter_cov_subdir = inter_subdirs["covariates_of_interest"]
        self.inter_tech_cov_subdir = inter_subdirs["technical_covariates"]
        inter_filenames = settings.get_setting("interaction_filenames")
        self.pvalue_filename = inter_filenames["pvalues"]
        self.zscore_filename = inter_filenames["zscores"]
        self.tvalue_filename = inter_filenames["tvalues"]

        self.eqtl_celltype = settings.get_setting("eqtl_celltype_datafile")
        self.celltypes = settings.get_setting("celltypes")
        self.colormap = settings.get_setting("colormap")
        self.cellmap_methods = settings.get_setting("cellmap_method_prefix_and_suffix")
        self.marker_genes = settings.get_setting("marker_genes_prefix")
        self.signif_cutoff = stats.norm.isf(settings.get_setting("significance_cutoff"))
        nrows = nrows
        if nrows == -1:
            nrows = None
        elif nrows <= 0:
            print("Unexpected argument for -n / --n_eqtls: '{}'".format(nrows))
            exit()
        self.nrows = nrows

        # Declare empty variables.
        self.eqtl_df = None
        self.geno_df = None
        self.alleles_df = None
        self.expr_df = None
        self.cov_df = None
        self.inter_cov_pvalue_df = None
        self.inter_tech_cov_pvalue_df = None
        self.inter_cov_zscore_df = None
        self.inter_tech_cov_zscore_df = None
        self.inter_cov_tvalue_df = None
        self.inter_tech_cov_tvalue_df = None
        self.eqtl_ct_df = None
        self.marker_df = None

    def load_all(self):
        print("Loading all dataframes for validation.")
        self.nrows = None
        self.get_celltypes()
        self.get_eqtl_df()
        self.get_geno_df()
        self.get_alleles_df()
        self.get_expr_df()
        self.get_cov_df()
        self.get_inter_cov_pvalue_df()
        self.get_inter_tech_cov_pvalue_df()
        self.get_inter_cov_zscore_df()
        self.get_inter_tech_cov_zscore_df()
        self.get_inter_cov_tvalue_df()
        self.get_inter_tech_cov_tvalue_df()
        self.get_eqtl_ct_df()
        self.get_marker_df()
        print("Validation finished.")
        exit()

    def get_celltypes(self):
        return self.celltypes

    def get_colormap(self):
        return self.colormap

    def get_cellmap_methods(self):
        return self.cellmap_methods

    def get_marker_genes(self):
        return self.marker_genes

    def get_significance_cutoff(self):
        return self.signif_cutoff

    def get_eqtl_df(self):
        if self.eqtl_df is None:
            self.eqtl_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                              self.eqtl_filename),
                                          header=0,
                                          index_col=False,
                                          nrows=self.nrows)
            self.eqtl_df.index = self.eqtl_df["SNPName"]
            self.eqtl_df.index.name = "-"
            self.validate()
        return self.eqtl_df

    def get_geno_df(self):
        if self.geno_df is None:
            self.geno_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                              self.geno_filename),
                                          header=0,
                                          index_col=0,
                                          nrows=self.nrows)
            self.validate()
        return self.geno_df

    def get_alleles_df(self):
        if self.alleles_df is None:
            self.alleles_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                                 self.alleles_filename),
                                             header=0,
                                             index_col=0,
                                             nrows=self.nrows)
            self.validate()
        return self.alleles_df

    def get_expr_df(self):
        if self.expr_df is None:
            self.expr_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                              self.expr_filename),
                                          header=0,
                                          index_col=0,
                                          nrows=self.nrows)
            self.validate()
        return self.expr_df

    def get_cov_df(self):
        if self.cov_df is None:
            self.cov_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                             self.cov_filename),
                                         header=0,
                                         index_col=0)
            self.validate()
        return self.cov_df

    def get_inter_cov_pvalue_df(self):
        if self.inter_cov_pvalue_df is None:
            self.inter_cov_pvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.pvalue_filename),
                header=0,
                index_col=0)
            self.validate()
        return self.inter_cov_pvalue_df

    def get_inter_tech_cov_pvalue_df(self):
        if self.inter_tech_cov_pvalue_df is None:
            self.inter_tech_cov_pvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.pvalue_filename),
                header=0,
                index_col=0)
            self.validate()
        return self.inter_tech_cov_pvalue_df

    def get_inter_cov_zscore_df(self):
        if self.inter_cov_zscore_df is None:
            self.inter_cov_zscore_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.zscore_filename),
                header=0,
                index_col=0)
            self.validate()
        return self.inter_cov_zscore_df

    def get_inter_tech_cov_zscore_df(self):
        if self.inter_tech_cov_zscore_df is None:
            self.inter_tech_cov_zscore_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.zscore_filename),
                header=0,
                index_col=0)
            self.validate()
        return self.inter_tech_cov_zscore_df

    def get_inter_cov_tvalue_df(self):
        if self.inter_cov_tvalue_df is None:
            self.inter_cov_tvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_cov_subdir,
                                    self.tvalue_filename),
                header=0,
                index_col=0)
            self.validate()
        return self.inter_cov_tvalue_df

    def get_inter_tech_cov_tvalue_df(self):
        if self.inter_tech_cov_tvalue_df is None:
            self.inter_tech_cov_tvalue_df = load_dataframe(
                inpath=os.path.join(self.inter_input_dir,
                                    self.inter_tech_cov_subdir,
                                    self.tvalue_filename),
                header=0,
                index_col=0)
            self.validate()
        return self.inter_tech_cov_tvalue_df

    def get_eqtl_ct_df(self):
        if self.eqtl_ct_df is None:
            self.eqtl_ct_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                                 self.eqtl_celltype),
                                             header=0,
                                             index_col=0)
            self.validate()
        return self.eqtl_ct_df

    def get_marker_df(self):
        if self.marker_df is None:
            self.marker_df = load_dataframe(inpath=os.path.join(self.input_dir,
                                                                self.markers_filename),
                                            header=0,
                                            index_col=False)
            self.validate()
        return self.marker_df

    def validate(self):
        dfs = [self.eqtl_df, self.geno_df, self.alleles_df, self.expr_df,
               self.eqtl_ct_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.index.identical(b.index):
                print("Order of eQTLs is not identical (1).")
                exit()

        dfs = [self.geno_df, self.expr_df, self.cov_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if a is not None and b is not None and \
                    not a.columns.identical(b.columns):
                print("Order of samples are not identical.")
                exit()

        print("\tValid.")
