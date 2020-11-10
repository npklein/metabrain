"""
File:         main.py
Created:      2020/10/13
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
from pathlib import Path
import itertools
import os

# Third party imports.

# Local application imports.
from local_settings import LocalSettings
from utilities import prepare_output_dir, check_file_exists, load_dataframe
from logger import Logger


class Main:
    def __init__(self, name, settings_file, clear_log):
        self.name = name
        self.settings_file = settings_file
        self.clear_log = clear_log

        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Get the arguments from the settings file.
        input_dir = settings.get_setting('input_dir')
        filenames = settings.get_setting('filenames')
        self.geno_file = os.path.join(input_dir, name, filenames["genotype"])
        self.expr_file = os.path.join(input_dir, name, filenames["expression"])
        self.cov_file = os.path.join(input_dir, name, filenames["covariates"])

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, name)
        prepare_output_dir(self.outdir)

        # Initialize logger.
        logger = Logger(outdir=self.outdir, clear_log=clear_log)
        self.log = logger.get_logger()

    def start(self):
        self.log.info("Starting program.")
        self.print_arguments()

        self.log.info("Loading data.")
        geno_df, expr_df, cov_df = self.load_data()

        self.log.info("Validating input data.")
        valid = self.validate(geno_df, expr_df, cov_df)
        if not valid:
            self.log.error("Input files do not match.")
            exit()
        else:
            self.log.info("Valid")

    def load_data(self):
        geno_df = None
        if check_file_exists(self.geno_file):
            geno_df = load_dataframe(inpath=self.geno_file,
                                     header=0,
                                     index_col=0,
                                     logger=self.log)
        else:
            self.log.error("Genotype file does not exist.")
            exit()

        expr_df = None
        if check_file_exists(self.expr_file):
            expr_df = load_dataframe(inpath=self.expr_file,
                                     header=0,
                                     index_col=0,
                                     logger=self.log)
        else:
            self.log.error("Expression file does not exist.")
            exit()

        cov_df = None
        if check_file_exists(self.cov_file):
            cov_df = load_dataframe(inpath=self.cov_file,
                                    header=0,
                                    index_col=0,
                                    logger=self.log)
        else:
            self.log.error("Covariate file does not exist.")
            exit()

        return geno_df, expr_df, cov_df

    @staticmethod
    def validate(geno_df, expr_df, cov_df):
        dfs = [geno_df, expr_df, cov_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if not a.columns.identical(b.columns):
                return False
        return True

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Name: {}".format(self.name))
        self.log.info("  > Settings file: {}".format(self.settings_file))
        self.log.info("  > Clear log: {}".format(self.clear_log))
        self.log.info("  > Genotype input file: {}".format(self.geno_file))
        self.log.info("  > Expression input file: {}".format(self.expr_file))
        self.log.info("  > Covariate input file: {}".format(self.cov_file))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("")