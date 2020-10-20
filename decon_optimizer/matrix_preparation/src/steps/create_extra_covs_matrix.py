"""
File:         create_extra_cov_matrices.py
Created:      2020/10/20
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
from functools import reduce
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe


class CreateExtraCovsMatrices:
    def __init__(self, settings, log, matrices, sample_dict, sample_order,
                 force, outdir):
        self.log = log
        self.matrices = matrices
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_extra_covs_matrices')
        prepare_output_dir(self.outdir)

    def start(self):
        self.log.info("Starting creating extra covariate file(s).")
        self.print_arguments()

        # Check if output file exist.
        for filepath in self.matrices:
            outpath = os.path.join(self.outdir, os.path.basename(filepath))
            if not check_file_exists(outpath) or self.force:
                self.log.info("Working on {}.".format(filepath))
                df = self.prepare_matrix(filepath)
                self.save(df, outpath)
            else:
                self.log.info("Skipping {}.".format(filepath))

    def prepare_matrix(self, filepath):
        self.log.info("\tLoading matrix.")
        df = load_dataframe(filepath, header=0, index_col=0, logger=self.log)

        self.log.info("\tPreprocessing.")
        df.columns = [self.sample_dict[x] if x in self.sample_dict else x for x in df.columns]

        self.log.info("\tChecking overlap.")
        overlap = [x for x in self.sample_order if x in df.columns]
        self.log.info("\t\t{}/{} [{:.2f}%] samples found".format(len(overlap), len(self.sample_order), (100/len(self.sample_order))*len(overlap)))

        self.log.info("\tSubsetting.")
        subset = df.loc[:, overlap]

        missing = set(self.sample_order) - set(subset.columns)
        if len(missing) > 0:
            self.log.info("\tCompleting data frame, adding {} missing sample columns.".format(len(missing)))
            missing_df = pd.DataFrame(np.nan, index=df.index, columns=missing)
            subset = subset.merge(missing_df, left_index=True, right_index=True)

        if list(subset.columns.to_list()) != self.sample_order:
            self.log.info("\tReordering columns.")
            subset = subset.loc[:, self.sample_order]

        return subset

    def save(self, df, outpath):
        print("\tSaving matrix.")
        save_dataframe(df=df, outpath=outpath,
                       index=True, header=True, logger=self.log)

    def clear_variables(self):
        self.matrices = None
        self.sample_dict = None
        self.sample_order = None
        self.force = None

    def print_arguments(self):
        self.log.info("Arguments:")
        self.log.info("  > Matrices: ")
        for i, file in enumerate(self.matrices):
            self.log.info("\t  > [{}] {}".format(i, file))
        self.log.info("  > Output directory: {}".format(self.outdir))
        self.log.info("  > Force: {}".format(self.force))
        self.log.info("")
