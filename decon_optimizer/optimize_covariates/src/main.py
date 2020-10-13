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


class Main:
    def __init__(self, name, settings_file):
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

    def start(self):
        print("Starting program.")
        self.print_arguments()

        print("Loading data.")
        geno_df, expr_df, cov_df = self.load_data()

        print("Validating input data.")
        self.validate(geno_df, expr_df, cov_df)

    def load_data(self):
        geno_df = None
        if check_file_exists(self.geno_file):
            geno_df = load_dataframe(inpath=self.geno_file,
                                     header=0,
                                     index_col=0)
        else:
            print("Genotype file does not exist.")
            exit()

        expr_df = None
        if check_file_exists(self.expr_file):
            expr_df = load_dataframe(inpath=self.expr_file,
                                     header=0,
                                     index_col=0)
        else:
            print("Expression file does not exist.")
            exit()

        cov_df = None
        if check_file_exists(self.cov_file):
            cov_df = load_dataframe(inpath=self.cov_file,
                                    header=0,
                                    index_col=0)
        else:
            print("Covariate file does not exist.")
            exit()

        return geno_df, expr_df, cov_df

    @staticmethod
    def validate(geno_df, expr_df, cov_df):
        dfs = [geno_df, expr_df, cov_df]
        for (a, b) in list(itertools.combinations(dfs, 2)):
            if not a.columns.identical(b.columns):
                print("Order of samples are not identical.")
                exit()
        print("\tValid")

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > Covariate input file: {}".format(self.cov_file))
        print("  > Output directory: {}".format(self.outdir))
        print("")