"""
File:         create_deconvolution_matrices.py
Created:      2020/04/07
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
import gzip
import os

# Third party imports.

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe


class CreateDeconvolutionMatrices:
    def __init__(self, settings, expr_file, expr_df, sample_dict, sample_order,
                 force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param expr_file: string, the expression data file.
        :param expr_df: DataFrame, the complete expression dataframe.
        :param sample_dict: dictionary, a dictionary for translating unmasked
                            sampels to the same format.
        :param sample_order: list, order of samples.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.celltype_profile_file = settings["celltype_profile_datafile"]
        self.translate_file = settings["translate_datafile"]
        self.marker_dict = settings["marker_dict"]
        self.expr_file = expr_file
        self.expr_df = expr_df
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_deconvolution_matrices')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.ct_profile_expr_outpath = os.path.join(self.outdir, "ct_profile_expr_table.txt.gz")
        self.markers_outpath = os.path.join(self.outdir, "marker_genes.txt.gz")

        # Create empty variable.
        self.celltype_profile = None

    def start(self):
        print("Starting creating matrices.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.markers_outpath) and \
                check_file_exists(self.ct_profile_expr_outpath) and \
                not self.force:
            print("Skipping step.")
            return

        # Load the complete expression file.
        if self.expr_df is None:
            # Load the expression matrix file.
            print("Loading expression matrix.")
            expr_df = load_dataframe(self.expr_file, header=0, index_col=0)
            expr_df = expr_df.rename(columns=self.sample_dict)
            self.expr_df = expr_df[self.sample_order]

        # Load the translate file.
        print("Loading translate matrix.")
        trans_df = load_dataframe(self.translate_file, header=0, index_col=None)
        trans_df = trans_df.loc[:, ["ArrayAddress", "Symbol"]]
        trans_df = trans_df.set_index("ArrayAddress")
        trans_df.index.name = "-"

        # Merge the two dataframes.
        complete_expr_df = self.expr_df.merge(trans_df, left_index=True, right_index=True)
        complete_expr_df = complete_expr_df.set_index("Symbol")
        complete_expr_df.index.name = "-"

        # Remove unneeded variables.
        self.expr_df = None
        del trans_df

        # Create the marker gene file.
        if not check_file_exists(self.markers_outpath) or self.force:
            if os.path.isfile(self.markers_outpath):
                print("Removing: {}".format(self.markers_outpath))
                os.remove(self.markers_outpath)

            print("Creating marker gene expression table.")
            marker_str_buffer = ["-" + "\t" + "\t".join(self.sample_order) + "\n"]
            for celltype, marker_genes in self.marker_dict.items():
                for marker_gene in marker_genes:
                    if marker_gene in complete_expr_df.index:
                        expression = complete_expr_df.loc[[marker_gene], :]
                        if (len(expression.index)) != 1:
                            print("\tMarker gene: {} gives 0 or >1 expression "
                                  "profiles.".format(marker_gene))
                            continue

                        marker_str = celltype + "_" + marker_gene + "\t" + \
                                     "\t".join(expression.iloc[0, :].astype(str).values) \
                                     + "\n"
                        marker_str_buffer.append(marker_str)
            self.write_buffer(self.markers_outpath, marker_str_buffer)

        # Create the marker gene file.
        if not check_file_exists(self.ct_profile_expr_outpath) or self.force:
            if os.path.isfile(self.ct_profile_expr_outpath):
                print("Removing: {}".format(self.ct_profile_expr_outpath))
                os.remove(self.ct_profile_expr_outpath)

            # Load the celltype profile file.
            print("Loading cell type profile matrix.")
            self.celltype_profile = load_dataframe(self.celltype_profile_file,
                                                   header=0, index_col=0)

            # Create the celltype profile file.
            print("Creating cell type profile expression table.")
            profile_str_buffer = ["-" + "\t" + "\t".join(self.sample_order) + "\n"]
            for marker_gene in self.celltype_profile.index:
                if marker_gene in complete_expr_df.index:
                    expression = complete_expr_df.loc[[marker_gene], :]
                    if (len(expression.index)) != 1:
                        print("\tMarker gene: {} gives 0 or >1 expression "
                              "profiles.".format(marker_gene))
                        continue

                    profile_str = marker_gene + "\t" + "\t".join(expression.iloc[0, :].astype(str).values) + "\n"
                    profile_str_buffer.append(profile_str)
            self.write_buffer(self.ct_profile_expr_outpath, profile_str_buffer)

    @staticmethod
    def write_buffer(filename, buffer):
        """
        Method for writing a list of strings to a gzipped file,

        :param filename: string, the output file path.
        :param buffer: list, the lines of strings to write.
        """
        # Write output files.
        if os.path.isfile(filename):
            mode = 'ab'
        else:
            mode = 'wb'

        with gzip.open(filename, mode) as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    def clear_variables(self):
        self.celltype_profile_file = None
        self.translate_file = None
        self.marker_dict = None
        self.expr_file = None
        self.expr_df = None
        self.sample_dict = None
        self.sample_order = None
        self.force = None

    def get_celltype_profile(self):
        return self.celltype_profile

    def get_celltype_profile_file(self):
        return self.celltype_profile_file

    def get_ct_profile_expr_outpath(self):
        return self.ct_profile_expr_outpath

    def get_markers_outpath(self):
        return self.markers_outpath

    def print_arguments(self):
        print("Arguments:")
        if self.expr_df is not None:
            print("  > Expression matrix: {}".format(self.expr_df.shape))
        else:
            print("  > Expression input file: {}".format(self.expr_file))
        print("  > Cell type profile input file: {}".format(self.celltype_profile_file))
        print("  > Cell type profile expression output path: {}".format(self.ct_profile_expr_outpath))
        print("  > Markers output path: {}".format(self.markers_outpath))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Force: {}".format(self.force))
        print("")
