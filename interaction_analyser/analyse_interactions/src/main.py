"""
File:         main.py
Created:      2020/03/13
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
import glob
import os

# Third party imports.

# Local application imports.
from src.utilities import get_project_root_dir, prepare_output_dir
from src.local_settings import LocalSettings
from src.utilities import get_extension, get_filename, check_file_exists


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, force, outdir):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param force: boolean, whether or not to force to redo each step.
        :param outdir: string, the name of the base output directory.
        """
        # Load the LocalSettings singelton class.
        settings = LocalSettings(settings_file)

        # Safe arguments.
        self.indir = settings.get_setting("input_dir")
        self.tech_covs = ' '.join(settings.get_setting("technical_covariates"))
        self.eqtl_ia = settings.get_setting("eQTLInteractionAnalyser")
        self.force = force

        # Prepare an output directory.
        self.ia_indir = os.path.join(get_project_root_dir(), outdir, 'input')
        prepare_output_dir(self.ia_indir)
        self.ia_outdir = os.path.join(get_project_root_dir(), outdir, 'output')
        prepare_output_dir(self.ia_outdir)

        # Construct filenames.
        self.eqtl_inpath = os.path.join(self.indir, 'eqtl_table.txt.gz')
        self.geno_inpath = os.path.join(self.indir, 'genotype_table.txt.gz')
        self.expr_inpath = os.path.join(self.indir, 'expression_table.txt.gz')
        self.cov_inpath = os.path.join(self.indir, 'covariates_table.txt.gz')

        # Cleanup output directory.
        self.cleanup(self.ia_outdir)
        if self.force:
            self.cleanup(self.ia_indir)

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
                    not check_file_exists(file3):
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
        print("Uncompressing / moving the eQTL file.\n")
        eqtl_inpath = self.uncompress_and_move(self.eqtl_inpath)

        # execute the program.
        print("\n### STEP3 ###\n")
        print("Executing eQTLInteractionAnalyser.\n")
        self.execute(eqtl_inpath)

        # remove the temporary files.
        print("\n### STEP4 ###\n")
        print("Removing uncompressed files.\n")
        self.cleanup(self.ia_indir, '.txt')

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

    @staticmethod
    def cleanup(dir, extension=None):
        if extension is None:
            extension = ""
        regex = "*" + extension

        for path in glob.glob(os.path.join(dir, regex)):
            if check_file_exists(path) and (
                    extension is None or get_extension(path) == extension):
                print("\r rm {}".format(path))
                os.remove(path)

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print(
            "  > eQTLInteractionAnalyser input directory: {}".format(self.ia_indir))
        print("  > eQTLInteractionAnalyser output directory: {}".format(
            self.ia_outdir))
        print("  > Technical covariates: {}".format(self.tech_covs))
        print("  > Interaction analyser: {}".format(self.eqtl_ia))
        print("  > Force: {}".format(self.force))
        print("")
