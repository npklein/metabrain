"""
File:         main.py
Created:      2020/11/16
Last Changed: 1010/12/21
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
from .cohort_object import Cohort
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
        decon_df = self.data.get_deco_fdr_df()
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
        del eqtl_decon_df

        self.log.info("Validating matrices.")
        self.validate(geno_df, expr_df, cf_df)

        self.log.info("Create cohort objects.")
        cohorts = {}
        sample_order = expr_df.columns
        for cohort, cohort_samples in cohort_sample_dict.items():
            mask = []
            for sample in sample_order:
                if sample not in cohort_samples:
                    mask.append(False)
                else:
                    mask.append(True)

            samples = set(geno_df.columns[mask])
            if len(samples) == 0:
                continue

            cohort_object = Cohort(cohort=cohort,
                                   samples=samples,
                                   geno_df=geno_df.loc[:, mask].copy(),
                                   expr_df=expr_df.loc[:, mask].copy(),
                                   cf_df=cf_df.loc[mask, :].copy(),
                                   log=self.log)

            #cohort_object.print_info()
            cohorts[cohort] = cohort_object
        del geno_df, expr_df, cf_df
        self.log.info("\tCreated {} cohort objects.".format(len(cohorts.keys())))

        self.log.info("Create cell type objects.")
        for cell_type in cell_types:
            if cell_type != "CellMapNNLS_Neuron":
                continue
            ct_mediated_eqtls = eqtl_signif_decon_df.loc[eqtl_signif_decon_df[cell_type] < self.alpha, :].copy()
            indices = ct_mediated_eqtls.index
            ct_mediated_eqtls.reset_index(drop=False, inplace=True)
            self.log.info("\tCell '{}' has {} cell type mediated eQTLs (FDR < {})".format(cell_type, ct_mediated_eqtls.shape[0], self.alpha))

            cell_type_object = CellType(cell_type=cell_type,
                                        eqtl_df=ct_mediated_eqtls,
                                        alleles_df=alleles_df.iloc[indices, :].copy(),
                                        indices=indices,
                                        cohorts=cohorts,
                                        sample_order=sample_order,
                                        log=self.log)

            cell_type_object.print_info()
            results = cell_type_object.test_all_cell_fractions(sample="HRA_01267")
            results.to_pickle(os.path.join(self.outdir, "results2.pkl"))
            exit()

        del alleles_df, eqtl_signif_decon_df

    def validate(self, geno_df, expr_df, ct_frac_df):
        if not ct_frac_df.index.equals(geno_df.columns):
            self.log.error("The genotype file columns do not match the cell "
                           "type fractions file.")
            exit()

        if not ct_frac_df.index.equals(expr_df.columns):
            self.log.error("The expressiom file columns do not match the cell "
                           "type fractions file.")
            exit()

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Alpha: {}".format(self.alpha))
        self.log.info("  > Extensions: {}".format(self.extensions))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("")




