"""
File:         main.py
Created:      2020/11/16
Last Changed: 1010/11/17
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
from __future__ import print_function
from pathlib import Path
import os

# Third party imports.

# Local application imports.
from .logger import Logger
from .data_object import Data
from .cell_type_object import CellType


class Main:
    def __init__(self, eqtl_path, genotype_path, alleles_path, expression_path,
                 cell_fractions_path, decon_path, sample_annotation_path,
                 sample_id, cohort_id, alpha, outdir, extensions):
        # Safe arguments.
        self.alpha = alpha
        self.extensions = extensions

        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, outdir)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Initialize logger.
        logger = Logger(outdir=self.outdir, clear_log=True)
        self.log = logger.get_logger()

        # Initialize data object.
        self.data = Data(eqtl_path=eqtl_path,
                         genotype_path=genotype_path,
                         alleles_path=alleles_path,
                         expression_path=expression_path,
                         cell_fractions_path=cell_fractions_path,
                         decon_path=decon_path,
                         sample_annotation_path=sample_annotation_path,
                         sample_id=sample_id,
                         cohort_id=cohort_id,
                         log=self.log)
        self.data.print_arguments()

    def start(self):
        self.log.info("Starting program.")
        self.print_arguments()

        self.log.info("Loading data")
        eqtl_df = self.data.get_eqtl_df()
        decon_df = self.data.get_deco_df()
        cell_types = [x for x in decon_df.columns if x not in ["SNPName", "ProbeName"]]
        self.log.info("\tCell types: {}".format(cell_types))

        eqtl_decon_df = eqtl_df.merge(decon_df,
                                      left_on=["SNPName", "ProbeName"],
                                      right_on=["SNPName", "ProbeName"])
        del eqtl_df, decon_df

        cf_df = self.data.get_frac_df()

        self.log.info("Filtering cell type mediated eQTLs")
        eqtl_signif_decon_df = eqtl_decon_df.loc[eqtl_decon_df[cell_types].min(axis=1) < self.alpha, :].copy()
        eqtl_signif_decon_df.reset_index(drop=True, inplace=True)

        self.log.info("Loading genotype / expression data of cell type mediated "
                      "eQTLs")
        skiprows = [x+1 for x in eqtl_decon_df.loc[eqtl_decon_df[cell_types].min(axis=1) >= self.alpha, :].index]
        geno_df = self.data.get_geno_df(skiprows=skiprows, nrows=max(eqtl_signif_decon_df.index))
        alleles_df = self.data.get_alle_df(skiprows=skiprows, nrows=max(eqtl_signif_decon_df.index))
        expr_df = self.data.get_expr_df(skiprows=skiprows, nrows=max(eqtl_signif_decon_df.index))
        cohort_sample_dict = self.data.get_cohort_to_sample_dict()

        self.log.info("Optimizing per cell type")
        for cell_type in cell_types:
            # Subset the cell type.
            ct_mediated_eqtls = eqtl_signif_decon_df.loc[eqtl_signif_decon_df[cell_type] < self.alpha, :].copy()
            ct_fractions = cf_df.loc[:, [cell_type in x for x in cf_df.columns]].copy()
            self.log.info("\tCell {} has {} cell type mediated eQTLs (FDR < {})".format(cell_type, ct_mediated_eqtls.shape[0], self.alpha))

            # Create the cell type object.
            self.log.info("\tCreating cell type object")
            cell_type = CellType(cell_type=cell_type,
                                 eqtl_df=ct_mediated_eqtls,
                                 geno_df=geno_df.iloc[ct_mediated_eqtls.index, :].copy(),
                                 alleles_df=alleles_df.iloc[ct_mediated_eqtls.index, :].copy(),
                                 expr_df=expr_df.iloc[ct_mediated_eqtls.index, :].copy(),
                                 ct_fractions=ct_fractions,
                                 cohort_sample_dict=cohort_sample_dict,
                                 log=self.log)
            results = cell_type.test_all_cell_fractions(sample="HRA_01267")
            results.to_pickle(os.path.join(self.outdir, "results.pkl"))
            exit()

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Alpha: {}".format(self.alpha))
        self.log.info("  > Extensions: {}".format(self.extensions))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("")




