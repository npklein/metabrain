"""
File:         main.py
Created:      2020/04/19
Last Changed: 2020/05/15
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
import scipy.stats as st

# Local application imports.
from general.utilities import prepare_output_dir
from general.df_utilities import load_dataframe, save_dataframe
from general.local_settings import LocalSettings


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Get the settings.
        self.eqtl_filepath = settings.get_setting("eqtl_datafile")
        self.interaction_filepath = settings.get_setting("interaction_datafile")
        self.marker_genes = settings.get_setting("marker_genes")
        self.marker_genes_min_true = settings.get_setting("marker_genes_min_true")
        self.decon_methods = settings.get_setting("deconvolution_methods")
        self.output_filepath = os.path.join(self.outdir, settings.get_setting("output_filename"))
        self.signif_cutoff = st.norm.isf(
            settings.get_setting("significance_cutoff"))

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Start assigning eQTLs to celltypes.")

        # Load the interaction matrix.
        inter_df = load_dataframe(inpath=self.interaction_filepath,
                                  header=0,
                                  index_col=0)

        # Load the eQTL data.
        eqtl_df = load_dataframe(inpath=self.eqtl_filepath,
                                 header=0,
                                 index_col=False,
                                 nrows=inter_df.shape[1])
        eqtl_df.index = eqtl_df["SNPName"]
        eqtl_df.index.name = None

        # Check if the files math up.
        for i in range(inter_df.shape[1]):
            if not inter_df.columns[i].startswith(eqtl_df.index[i]):
                print("Input files do not match.")
                exit()

        # Loop over the marker genes.
        marker_data = inter_df.loc[inter_df.index.str.startswith(self.marker_genes), :]
        celltypes = set([x.split("_")[1] for x in marker_data.index])
        for celltype in celltypes:
            celltype_data = marker_data.loc[marker_data.index.str.startswith(self.marker_genes + celltype), :]
            signif_inter = (celltype_data > abs(self.signif_cutoff)).astype(int)
            eqtl_df.loc[:, self.marker_genes + celltype] = (signif_inter.sum() >= self.marker_genes_min_true).astype(int)

        # Loop over each deconvolution method.
        all_methods = self.decon_methods.copy()
        all_methods.append(self.marker_genes)
        for prefix in all_methods:
            method_data = inter_df.loc[inter_df.index.str.startswith(prefix), :]

            for (index, z_scores) in method_data.iterrows():
                eqtl_df.loc[:, index] = (z_scores > abs(self.signif_cutoff)).astype(int)

        # Save the new dataframe.
        save_dataframe(eqtl_df, self.output_filepath, header=True, index=False)

    def print_arguments(self):
        """
        Method for printing the variables of the class.
        """
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > EQTL data file: {}".format(self.eqtl_filepath))
        print("  > Interaction data file: {}".format(self.interaction_filepath))
        print("  > Marker genes: {}".format(self.marker_genes))
        print("  > Marker genes min. True: {}".format(self.marker_genes_min_true))
        print("  > Deconvolution Methods: {}".format(self.decon_methods))
        print("  > Output file: {}".format(self.output_filepath))
        print("")
