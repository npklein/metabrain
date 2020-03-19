"""
File:         main.py
Created:      2020/03/19
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
import pickle
import glob
import os
import re

# Third party imports.
import numpy as np

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.local_settings import LocalSettings
from general.utilities import get_leaf_dir, get_basename
from general.df_utilities import load_dataframe, save_dataframe


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, groups, force):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param groups: list, the names of groups to analyse.
        :param force: boolean, whether or not to force to redo each step.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.data_indir = settings.get_setting("data_dir")
        self.g_data_indir = settings.get_setting("groups_data_dir")
        self.g_inter_indir = settings.get_setting("inter_groups_dir")
        self.groups = groups
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Find which groups are in both input directories.
        self.group_ids = self.filter_groups()

        # Prepare filenames.
        filenames = settings.get_setting("filenames")
        self.obj_filename = filenames["object"]
        self.eqtl_filename = filenames["eqtl"]
        self.geno_filename = filenames["genotype"]
        self.expr_filename = filenames["expression"]
        self.cov_filename = filenames["covariate"]
        self.inter_regex = filenames["interaction"]

    def filter_groups(self):
        g_in_data_indir = glob.glob(os.path.join(self.g_data_indir, 'group_*'))
        g_data_ids = set([get_leaf_dir(x) for x in g_in_data_indir])

        g_in_inter_indir = glob.glob(os.path.join(self.g_inter_indir,'group_*'))
        g_inter_ids = set([get_leaf_dir(x) for x in g_in_inter_indir])

        return g_data_ids.intersection(g_inter_ids)

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting combining groups.")
        self.print_arguments()

        # Define output filenames.
        eqtl_outpath = os.path.join(self.outdir, "interaction_table.txt.gz")
        geno_outpath = os.path.join(self.outdir, "genotype_table.txt.gz")
        alleles_outpath = os.path.join(self.outdir, "alleles_table.txt.gz")
        expr_outpath = os.path.join(self.outdir, "expression_table.txt.gz")
        cov_outpath = os.path.join(self.outdir, "covariate_table.txt.gz")
        inter_outpath = os.path.join(self.outdir, "interaction_table.txt.gz")

        snp_indices = np.array([], dtype=np.int16)
        sample_indices = np.array([], dtype=np.int16)
        inter_df = None

        print("Loading group data.")
        for i, group_id in enumerate(self.group_ids):
            print("\tWorking on: {:10s} [{}/{} "
                  "{:.2f}%]".format(group_id, i + 1, len(self.group_ids),
                                    (100 / len(self.group_ids)) * i + 1))

            # Define the directory names.
            data_indir = os.path.join(self.g_data_indir, group_id)
            inter_indir = os.path.join(self.g_inter_indir, group_id, 'output')

            # Load the group object.
            obj_inpath = os.path.join(data_indir, self.obj_filename + ".pkl")
            with open(obj_inpath, "rb") as f:
                group_object = pickle.load(f)

            # Safe the indices.
            snp_indices = np.append(snp_indices, group_object.get_snp_indices())
            sample_indices = np.append(snp_indices, group_object.get_sample_indices())

            if not check_file_exists(inter_outpath) or self.force:
                # Search for the interaction filename.
                inter_inpath = None
                for path in glob.glob(os.path.join(inter_indir, "*")):
                    if re.match(self.inter_regex, get_basename(path)):
                        inter_inpath = path
                        break
                if inter_inpath is None:
                    print("Interaction matrix not found.")
                    exit()

                # Load the interaction file.
                group_inter_df = load_dataframe(inpath=inter_inpath,
                                                header=0,
                                                index_col=0)

                # Merge them.
                if inter_df is None:
                    inter_df = group_inter_df
                else:
                    inter_df = inter_df.merge(group_inter_df,
                                              left_index=True,
                                              right_index=True)

        if not check_file_exists(inter_outpath) or self.force:
            print("Constructing complete interaction matrix.")

            # Set the indices.
            inter_df = inter_df.T
            inter_df["index"] = snp_indices
            inter_df.sort_values(by=['index'], inplace=True)
            inter_df.drop(["index"], axis=1, inplace=True)
            inter_df = inter_df.T

            # Safe.
            save_dataframe(df=inter_df, outpath=inter_outpath,
                           index=True, header=True)

    def print_arguments(self):
        print("Arguments:")
        print("  > Data input directory: {}".format(self.data_indir))
        print("  > Group data input directory: {}".format(self.g_data_indir))
        print("  > Group interaction input directory: {}".format(self.g_inter_indir))
        print("  > Group ids: {}".format(self.group_ids))
        print("  > Force: {}".format(self.force))
        print("")
