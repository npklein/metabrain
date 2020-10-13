"""
File:         main.py
Created:      2020/10/08
Last Changed: 2020/10/13
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
from local_settings import LocalSettings
from utilities import prepare_output_dir
from logger import Logger
from .steps.combine_gte_files import CombineGTEFiles
from .steps.combine_eqtlprobes import CombineEQTLProbes
from .steps.create_cohort_matrix import CreateCohortMatrix
from .steps.create_matrices import CreateMatrices
from .steps.correct_cohort_effects import CorrectCohortEffects
from .steps.perform_deconvolution import PerformDeconvolution
from .steps.filter_technical_covariates import FilterTechnicalCovariates
from .steps.create_cov_matrix import CreateCovMatrix


class Main:
    def __init__(self, name, settings_file, force_steps, clear_log):
        self.name = name
        self.settings_file = settings_file
        self.force_steps = force_steps
        self.clear_log = clear_log

        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.force_dict = self.create_force_dict(force_steps)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, name)
        prepare_output_dir(self.outdir)

        # Initialize logger.
        logger = Logger(outdir=self.outdir, clear_log=clear_log)
        self.log = logger.get_logger()

    @staticmethod
    def create_force_dict(force_steps):
        step_dependencies = {'combine_gte_files': ['create_cohort_matrix', 'create_matrices', 'correct_cohort_effects', 'perform_deconvolution', 'create_cov_matrix'],
                             'combine_eqtlprobes': ['create_cohort_matrix', 'create_matrices', 'correct_cohort_effects', 'perform_deconvolution', 'create_cov_matrix'],
                             'create_cohort_matrix': ['create_matrices', 'correct_cohort_effects', 'perform_deconvolution', 'create_cov_matrix'],
                             'create_matrices': ['correct_cohort_effects', 'perform_deconvolution', 'create_cov_matrix'],
                             'correct_cohort_effects': ['perform_deconvolution', 'create_cov_matrix'],
                             'perform_deconvolution': ['create_cov_matrix'],
                             'filter_technical_covariates': ['create_cov_matrix'],
                             'create_cov_matrix': []}
        force_dict = {step: False for step in step_dependencies.keys()}

        if force_steps is None or len(force_steps) == 0:
            return force_dict

        if force_steps == ['all']:
            for key in force_dict.keys():
                force_dict[key] = True
        else:
            for step in force_steps:
                dependencies = step_dependencies[step]
                dependencies.append(step)
                for substep in dependencies:
                    force_dict[substep] = True

        return force_dict

    def start(self):
        self.log.info("Starting program.")
        self.print_arguments()

        self.log.info("### STEP1 ###")
        self.log.info("")
        # Step 1. Combine GTE files.
        cgtef = CombineGTEFiles(
            settings=self.settings.get_setting('combine_gte_files'),
            log=self.log,
            force=self.force_dict['combine_gte_files'],
            outdir=self.outdir)
        cgtef.start()
        cgtef.clear_variables()
        self.log.info("")

        # Step2. Combine eQTL probes files.
        self.log.info("### STEP2 ###")
        self.log.info("")
        cepf = CombineEQTLProbes(
            settings=self.settings.get_setting('combine_eqtlprobes'),
            log=self.log,
            force=self.force_dict['combine_eqtlprobes'],
            outdir=self.outdir)
        cepf.start()
        cepf.clear_variables()
        self.log.info("")

        # Step3. Create the cohort matrix.
        self.log.info("### STEP3 ###")
        self.log.info("")
        ccm = CreateCohortMatrix(
            settings=self.settings.get_setting('create_cohort_matrix'),
            log=self.log,
            sample_dict=cgtef.get_sample_dict(reverse=True),
            sample_order=cgtef.get_sample_order(),
            force=self.force_dict['create_cohort_matrix'],
            outdir=self.outdir)
        ccm.start()
        ccm.clear_variables()
        self.log.info("")

        # Step4. Create the ordered matrices.
        self.log.info("### STEP4 ###")
        self.log.info("")
        cm = CreateMatrices(
            settings=self.settings.get_setting('create_matrices'),
            log=self.log,
            sample_dict=cgtef.get_sample_dict(),
            sample_order=cgtef.get_sample_order(),
            eqtl_file=cepf.get_eqtl_file(),
            eqtl_df=cepf.get_eqtl_df(),
            force=self.force_dict['create_matrices'],
            outdir=self.outdir)
        cm.start()
        cm.clear_variables()
        self.log.info("")

        # Step5. Correct the gene expression for cohort effects.
        self.log.info("### STEP5 ###")
        self.log.info("")
        cce = CorrectCohortEffects(
            settings=self.settings.get_setting('correct_cohort_effects'),
            log=self.log,
            cohort_file=ccm.get_cohort_file(),
            cohort_df=ccm.get_cohort_df(),
            expr_file=cm.get_expr_file(),
            expr_df=cm.get_expr_df(),
            sign_expr_file=cm.get_sign_expr_file(),
            sign_expr_df=cm.get_sign_expr_df(),
            force=self.force_dict['correct_cohort_effects'],
            outdir=self.outdir)
        cce.start()
        cce.clear_variables()
        self.log.info("")

        # Step6. Perform NNLS deconvolution.
        self.log.info("### STEP6 ###")
        self.log.info("")
        pd = PerformDeconvolution(
            settings=self.settings.get_setting('perform_deconvolution'),
            log=self.log,
            sign_file=cm.get_sign_file(),
            sign_df=cm.get_sign_df(),
            sign_expr_file=cce.get_sign_expr_cc_file(),
            sign_expr_df=cce.get_sign_expr_cc_df(),
            force=self.force_dict['perform_deconvolution'],
            outdir=self.outdir)
        pd.start()
        pd.clear_variables()
        self.log.info("")

        # Step7. Filter technical covariates.
        self.log.info("### STEP7 ###")
        self.log.info("")
        ftc = FilterTechnicalCovariates(
            settings=self.settings.get_setting('filter_technical_covariates'),
            log=self.log,
            force=self.force_dict['filter_technical_covariates'],
            sample_dict=cgtef.get_sample_dict(),
            sample_order=cgtef.get_sample_order(),
            outdir=self.outdir)
        ftc.start()
        ftc.clear_variables()
        self.log.info("")

        # Step8. Create the covariance matrix.
        self.log.info("### STEP8 ###")
        self.log.info("")
        ccm = CreateCovMatrix(
            settings=self.settings.get_setting('create_cov_matrix'),
            log=self.log,
            tech_covs_file=ftc.get_tech_covs_file(),
            tech_covs_df=ftc.get_tech_covs_df(),
            cohort_file=ccm.get_cohort_file(),
            cohort_df=ccm.get_cohort_df(),
            decon_file=pd.get_decon_file(),
            decon_df=pd.get_decon_df(),
            sample_dict=cgtef.get_sample_dict(),
            sample_order=cgtef.get_sample_order(),
            force=self.force_dict['create_cov_matrix'],
            outdir=self.outdir)
        ccm.start()
        ccm.clear_variables()
        self.log.info("")

        # End.
        self.log.info("")
        self.log.info("Program completed.")

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Name: {}".format(self.name))
        self.log.info("  > Settings file: {}".format(self.settings_file))
        self.log.info("  > Force steps: {}".format(self.force_steps))
        self.log.info("  > Clear log: {}".format(self.clear_log))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("")
