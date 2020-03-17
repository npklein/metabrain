"""
File:         main.py
Created:      2020/03/13
Last Changed: 2020/03/17
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
import glob
import os

# Third party imports.
import pandas as pd
import scipy.stats as st

# Local application imports.
from src.utilities import get_project_root_dir, prepare_output_dir
from src.local_settings import LocalSettings
from src.utilities import get_extension, get_filename, check_file_exists
from src.df_utilities import load_dataframe, save_dataframe


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, force):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param force: boolean, whether or not to force to redo each step.
        """
        # Load the LocalSettings singelton class.
        settings = LocalSettings(settings_file)

        # Safe arguments.
        self.indir = settings.get_setting("input_dir")
        self.tech_covs = ' '.join(settings.get_setting("technical_covariates"))
        self.eqtl_ia = settings.get_setting("eQTLInteractionAnalyser")
        self.celltypes = settings.get_setting("celltypes")
        self.force = force

        # Prepare an output directory.
        outdir = settings.get_setting('output_dir')
        self.ia_indir = os.path.join(get_project_root_dir(), outdir, 'input')
        prepare_output_dir(self.ia_indir)
        self.ia_outdir = os.path.join(get_project_root_dir(), outdir, 'output')
        prepare_output_dir(self.ia_outdir)

        # Construct filenames.
        self.eqtl_inpath = os.path.join(self.indir, 'eqtl_table.txt.gz')
        self.geno_inpath = os.path.join(self.indir, 'genotype_table.txt.gz')
        self.expr_inpath = os.path.join(self.indir, 'expression_table.txt.gz')
        self.cov_inpath = os.path.join(self.indir, 'covariates_table.txt.gz')

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        # EQTLInteractioAnalyser expected input.
        print("\n### STEP1 ###\n")
        expected_input = ["Genotypes", "Expression", "Covariates"]
        for filename in expected_input:
            file1 = os.path.join(self.ia_indir, filename + ".binary.dat")
            file2 = os.path.join(self.ia_indir, filename + ".binary.rows.txt")
            file3 = os.path.join(self.ia_indir,
                                 filename + ".binary.columns.txt")

            if not check_file_exists(file1) or \
                    not check_file_exists(file2) or \
                    not check_file_exists(file3) or \
                    self.force:
                # Uncompressing the input files.
                print("Uncompressing / moving input files.\n")
                geno_inpath = self.uncompress_and_move(self.geno_inpath)
                expr_inpath = self.uncompress_and_move(self.expr_inpath)
                cov_inpath = self.uncompress_and_move(self.cov_inpath)

                # Convert to binary.
                print("Converting files to binary format.\n")
                self.convert_to_binary(geno_inpath, 'Genotypes')
                self.convert_to_binary(expr_inpath, 'Expression')
                self.convert_to_binary(cov_inpath, 'Covariates')
            else:
                print("Skipping {} preparation.".format(filename))

        print("\n### STEP2 ###\n")
        eqtl_inpath = os.path.join(self.ia_indir,
                                   get_filename(self.eqtl_inpath) + ".txt")
        if not check_file_exists(eqtl_inpath) or self.force:
            print("Uncompressing / moving the eQTL file.\n")
            eqtl_inpath = self.uncompress_and_move(self.eqtl_inpath)
        else:
            print("Skipping step.")

        # execute the program.
        print("\n### STEP3 ###\n")
        inter_files = glob.glob(
            os.path.join(self.ia_outdir, "InteractionZScoresMatrix*"))
        if len(inter_files) != 1 or self.force:
            print("Executing the eQTLInteractionAnalyser.\n")
            self.execute(eqtl_inpath)
        else:
            print("Skipping step.")

        # perform eQTL interaction clustering on marker genes.
        print("\n### STEP4 ###\n")
        inter_inpath = glob.glob(os.path.join(self.ia_outdir, "InteractionZScoresMatrix*"))[0]
        self.eqtl_celltype_mediated(inter_inpath)

    def uncompress_and_move(self, inpath):
        outpath = os.path.join(self.ia_indir, get_filename(inpath) + ".txt")

        if get_extension(inpath) == '.txt.gz':
            command = 'gunzip -c {} > {}'.format(inpath, outpath)
            print("\t{}".format(command))
            os.system(command)
        else:
            command = 'cp {} {}'.format(inpath, outpath)
            print("\t{}".format(command))
            os.system(command)

        return outpath

    def convert_to_binary(self, inpath, out_filename):
        outpath = os.path.join(self.ia_indir, out_filename + ".binary")

        command = 'java -jar {} --convertMatrix -i {} -o {}'.format(
            self.eqtl_ia, inpath, outpath)
        print("\t{}".format(command))
        os.system(command)

    def execute(self, eqtl_inpath):
        command = 'java -jar {} -input {} -output {} -eqtls {} ' \
                  '-maxcov 1 ' \
                  '-noNormalization ' \
                  '-nnoCovNormalization ' \
                  '-cov {}'.format(self.eqtl_ia, self.ia_indir, self.ia_outdir,
                                   eqtl_inpath, self.tech_covs)
        print("\t{}".format(command))
        os.system(command)

    def eqtl_celltype_mediated(self, inter_inpath):
        # Loading dataframes.
        print("Loading interaction dataframe.")
        inter_df = load_dataframe(inpath=inter_inpath,
                                  header=0,
                                  index_col=0)
        inter_df = inter_df.T

        print("Loading eQTL dataframe.")
        eqtl_df = load_dataframe(inpath=self.eqtl_inpath,
                                 header=0,
                                 index_col=1)
        eqtl_df = eqtl_df[["ProbeName", "HGNCName"]]

        # Calculate the z-score cutoff.
        z_score_cutoff = st.norm.ppf(
            0.05 / (inter_df.shape[0] * inter_df.shape[1]))
        gini_cutoff = 0.75

        # Subset on the marker genes.
        marker_cols = []
        for colname in inter_df.columns:
            if ("_" in colname) and (colname.split("_")[0] in self.celltypes):
                marker_cols.append(colname)

        marker_df = inter_df.loc[:, marker_cols]
        del inter_df

        # Create a gini dataframe grouped by celltype.
        gini_df = marker_df.copy()
        gini_df = gini_df.abs()
        zscore_mask = list(gini_df.max(axis=1) >= abs(z_score_cutoff))
        gini_df.columns = [x.split("_")[0] for x in gini_df.columns]
        gini_df = gini_df.T.groupby(gini_df.columns).sum().T

        # Calculate the gini impurity.
        gini_values = gini_df.div(gini_df.sum(axis=1), axis=0).pow(2)
        marker_df["gini_impurity"] = 1 - gini_values.sum(axis=1)
        marker_df["eqtl_celltype"] = gini_values.idxmax(axis=1)

        # Subset the marker df on gini impurity.
        gini_mask = list(marker_df["gini_impurity"] <= gini_cutoff)
        marker_df = marker_df.loc[zscore_mask and gini_mask, :]
        marker_df.index.name = "-"
        marker_df.reset_index(inplace=True)

        # Subset the eQTL dataframe.
        eqtl_df = eqtl_df.loc[zscore_mask and gini_mask, :]
        eqtl_df.reset_index(inplace=True)

        # Merge them together.
        merged_df = pd.concat([marker_df, eqtl_df], axis=1)
        merged_df = merged_df.sort_values(by=['eqtl_celltype', 'gini_impurity'])

        # Save the dataframe.
        save_dataframe(df=merged_df,
                       outpath=os.path.join(self.ia_outdir,
                                            "marker_table.txt.gz"),
                       header=True,
                       index=False)

        # Save celltype eqtl's HGNC names.
        print("Writing celltype mediated eQTL files.")
        for celltype in marker_df['eqtl_celltype'].unique():
            subset = merged_df.loc[merged_df['eqtl_celltype'] == celltype, :]
            print("\tCelltype: {}\t{} genes".format(celltype, len(subset.index)))
            if len(subset.index) > 0:
                genes = ', '.join(subset['HGNCName'].to_list())
                outfile = open(os.path.join(self.ia_outdir,
                                            '{}.txt'.format(celltype)), "w")
                outfile.write(genes)
                outfile.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print(
            "  > eQTLInteractionAnalyser input directory: {}".format(
                self.ia_indir))
        print("  > eQTLInteractionAnalyser output directory: {}".format(
            self.ia_outdir))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > Interaction analyser: {}".format(self.eqtl_ia))
        print("  > Force: {}".format(self.force))
        print("")
