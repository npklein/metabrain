"""
File:         main.py
Created:      2020/03/12
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
from __future__ import print_function
import os

# Third party imports.

# Local application imports.
from src.utilities import get_project_root_dir, prepare_output_dir
from src.df_utilities import load_dataframe
from src.local_settings import LocalSettings
from src.steps.combine_gte_files import CombineGTEFiles
from src.steps.combine_eqtlprobes import CombineEQTLProbes
from src.steps.create_matrices import CreateMatrices
from src.steps.create_cov_matrix import CreateCovMatrix
from src.steps.mask_matrices import MaskMatrices
from src.steps.create_groups import CreateGroups


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, force_steps):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param force_steps: list, the names of the steps to force to redo.
        """
        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(settings_file)

        # Safe arguments.
        self.force_dict = self.create_force_dict(force_steps)

        # Prepare an output directory.
        self.outdir = os.path.join(get_project_root_dir(), 'output')
        prepare_output_dir(self.outdir)

    @staticmethod
    def create_force_dict(force_steps):
        force_dict = {'combine_gte_files': False, 'combine_eqtlprobes': False,
                      'create_matrices': False, 'create_cov_matrix': False,
                      'mask_matrices': False, 'create_groups': False}
        if force_steps is None or len(force_steps) == 0:
            return force_dict

        if force_steps == ['all']:
            for key in force_dict.keys():
                force_dict[key] = True
        else:
            for step in force_steps:
                if step in force_dict.keys():
                    force_dict[step] = True

        return force_dict

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting program.")
        print("\n### STEP1 ###\n")
        # Step 1. Combine GTE files.
        cgtef = CombineGTEFiles(
            settings=self.settings.get_setting('combine_gte_files'),
            force=self.force_dict['combine_gte_files'],
            outdir=self.outdir)
        cgtef.start()

        # Step2. Combine eQTL probes files.
        print("\n### STEP2 ###\n")
        cepf = CombineEQTLProbes(
            settings=self.settings.get_setting('combine_eqtlprobes'),
            force=self.force_dict['combine_eqtlprobes'],
            outdir=self.outdir)
        cepf.start()

        # Step3. Create the ordered unmasked matrices.
        print("\n### STEP3 ###\n")
        cm = CreateMatrices(
            settings=self.settings.get_setting('create_matrices'),
            gte_df=cgtef.get_gte(),
            sample_dict=cgtef.get_sample_dict(),
            sample_order=cgtef.get_sample_order(),
            eqtl_df=cepf.get_eqtlprobes(),
            force=self.force_dict['create_matrices'],
            outdir=self.outdir)
        cm.start()

        # Step4. Create the covariance matrix.
        print("\n### STEP4 ###\n")
        ccm = CreateCovMatrix(
            settings=self.settings.get_setting('create_matrices'),
            marker_file=cm.get_markers_outpath(),
            force=self.force_dict['create_matrices'],
            outdir=self.outdir)
        ccm.start()

        # Load the complete dataframes.
        print("\n### LOADING DATAFRAMES ###\n")
        geno_df = load_dataframe(cm.get_geno_outpath(),
                                 header=True,
                                 index_col=True)
        alleles_df = load_dataframe(cm.get_alleles_outpath(),
                                    header=True,
                                    index_col=True)
        expr_df = load_dataframe(cm.get_expr_outpath(),
                                 header=True,
                                 index_col=True)
        cov_df = ccm.get_covariates()

        # Step5. Create the masked matrices.
        print("\n### STEP5 ###\n")
        cmm = MaskMatrices(
            settings=self.settings.get_setting('mask_matrices'),
            geno_df=geno_df.copy(),
            alleles_df=alleles_df.copy(),
            expr_df=expr_df.copy(),
            cov_df=cov_df.copy(),
            sample_dict_masked=cgtef.get_sample_dict(),
            eqtl_dict_masked=cepf.get_eqtl_dict_masked(),
            force=self.force_dict['mask_matrices'],
            outdir=self.outdir)
        cmm.start()

        # Step 6. Create the group matrices.
        print("\n### STEP6 ###\n")
        cg = CreateGroups(
            settings=self.settings.get_setting('create_groups'),
            geno_df=geno_df.copy(),
            alleles_df=alleles_df.copy(),
            expr_df=expr_df.copy(),
            cov_df=cov_df.copy(),
            groups_file=cm.get_group_outpath(),
            force=self.force_dict['create_groups'],
            outdir=self.outdir)
        cg.start()

