"""
File:         main.py
Created:      2020/03/13
Last Changed: 2020/03/19
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
import subprocess
import glob
import os
import re

# Third party imports.

# Local application imports.
from general.local_settings import LocalSettings
from general.utilities import get_leaf_dir, check_file_exists, get_basename, \
    prepare_output_dir


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, groups, force, verbose):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param groups: list, the names of groups to analyse.
        :param force: boolean, whether or not to force to redo each step.
        :param verbose: boolean, whether or not to print each step.
        """
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Safe arguments.
        self.indir = settings.get_setting("input_dir")
        self.tech_covs = settings.get_setting("technical_covariates")
        self.eqtl_ia = settings.get_setting("eQTLInteractionAnalyser")
        self.groups = groups
        self.force = force
        self.verbose = verbose

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   settings.get_setting("output_dir"))
        prepare_output_dir(self.outdir)

        # Find which groups are in the input directory.
        groups_in_indir = glob.glob(os.path.join(self.indir, 'group_*'))
        self.group_indirs = self.filter_groups(groups_in_indir)

        # Prepare filenames.
        filenames = settings.get_setting("filenames")
        self.eqtl_filename = filenames["eqtl"]
        self.geno_filename = filenames["genotype"]
        self.expr_filename = filenames["expression"]
        self.cov_filename = filenames["covariate"]

    def filter_groups(self, groups_in_indir):
        # Get the groups we are interested in.
        group_indirs = []
        if "all" in self.groups:
            group_indirs = groups_in_indir
        else:
            for group_id in self.groups:
                if re.match("^group_[0-9]+$", group_id):
                    group_indir = os.path.join(self.indir, group_id)
                    if group_indir in groups_in_indir:
                        group_indirs.append(group_indir)
                    else:
                        print("Unexpected group_id: {}.".format(group_id))

        return group_indirs

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting interaction analyser.")
        self.print_arguments()

        # Loop over the groups.
        print("Performing interaction analyses.")
        for i, group_indir in enumerate(self.group_indirs):
            group_id = get_leaf_dir(group_indir)
            print("\tWorking on: {:10s} [{}/{} "
                  "{:.2f}%]".format(group_id, i + 1, len(self.group_indirs),
                                    (100 / len(self.group_indirs)) * i + 1))

            # Prepare the input and output directories.
            group_outdir = os.path.join(self.outdir, group_id)
            ia_indir = os.path.join(group_outdir, 'input')
            ia_outdir = os.path.join(group_outdir, 'output')
            for outdir in [group_outdir, ia_indir, ia_outdir]:
                prepare_output_dir(outdir)

            # Prepare the EQTLInteractioAnalyser expected input.
            self.print_string("\n### STEP1 ###\n")
            expected_input = ["Genotypes", "Expression", "Covariates"]
            filenames = [self.geno_filename, self.expr_filename,
                         self.cov_filename]
            for exp_ia_infile, filename in zip(expected_input, filenames):
                # Check if the files alreadt exist.
                file1 = os.path.join(ia_indir, exp_ia_infile + ".binary.dat")
                file2 = os.path.join(ia_indir,
                                     exp_ia_infile + ".binary.rows.txt")
                file3 = os.path.join(ia_indir,
                                     exp_ia_infile + ".binary.columns.txt")

                if not check_file_exists(file1) or \
                        not check_file_exists(file2) or \
                        not check_file_exists(file3) or \
                        self.force:
                    self.print_string("\nPreparing {}.".format(filename))

                    # Define the filenames.
                    compr_file = os.path.join(self.indir, group_id,
                                              filename + '.txt.gz')
                    copy_file = os.path.join(ia_indir, filename + '.txt.gz')
                    uncompr_file = os.path.join(ia_indir, filename + '.txt')
                    bin_file = os.path.join(ia_indir, exp_ia_infile + ".binary")

                    # Copy and decompressed the file.
                    self.print_string("\nCopying the input files.")
                    self.copy_file(compr_file, copy_file)
                    self.print_string("\nDecompressing the input files.")
                    self.decompress(copy_file)

                    # Convert to binary.
                    self.print_string("\nConverting files to binary format.")
                    self.convert_to_binary(uncompr_file, bin_file)

                    # Remove the uncompressed file.
                    self.print_string("\nRemoving uncompressed files.")
                    if check_file_exists(uncompr_file):
                        self.print_string("\tos.remove({})".format(uncompr_file))
                        os.remove(uncompr_file)
                else:
                    self.print_string("Skipping {} preparation.".format(filename))

            # prepare the eQTL file.
            self.print_string("\n### STEP2 ###\n")
            eqtl_file = os.path.join(ia_indir, self.eqtl_filename + '.txt')
            if not check_file_exists(eqtl_file) or self.force:
                self.print_string("\nPreparing eQTL file.")
                # Define the filenames.
                compr_file = os.path.join(self.indir, group_id,
                                          self.eqtl_filename + '.txt.gz')
                copy_file = os.path.join(ia_indir, self.eqtl_filename + '.txt.gz')

                # Copy and decompressed the file.
                self.print_string("\nCopying the input files.")
                self.copy_file(compr_file, copy_file)
                self.print_string("\nDecompressing the input files.")
                self.decompress(copy_file)
            else:
                self.print_string("Skipping eqtl preparation.")

            # execute the program.
            self.print_string("\n### STEP3 ###\n")
            found = False
            if not self.force:
                regex = "^InteractionZScoresMatrix-[0-9]+Covariates.txt$"
                for path in glob.glob(os.path.join(ia_outdir, "*")):
                    if re.match(regex, get_basename(path)):
                        found = True
                        break

            if self.force or not found:
                self.print_string("Executing the eQTLInteractionAnalyser.")
                self.execute(ia_indir, ia_outdir, eqtl_file)

    def print_string(self, string):
        if self.verbose:
            print(string)

    def copy_file(self, inpath, outpath):
        command = ['cp', inpath, outpath]
        self.execute_command(command)

    def decompress(self, inpath):
        command = ['gunzip', inpath, '-f']
        self.execute_command(command)

    def convert_to_binary(self, inpath, outpath):
        command = ['java', '-jar', self.eqtl_ia,
                   '--convertMatrix',
                   '-i', inpath, '-o', outpath]
        self.execute_command(command)

    def execute(self, ia_indir, ia_outdir, eqtl_inpath):
        command = ['java', '-jar', self.eqtl_ia,
                   '-input', ia_indir, '-output', ia_outdir,
                   '-eqtls', eqtl_inpath,
                   '--maxcov', '1',
                   '-noNormalization',
                   '-noCovNormalization',
                   '-cov', *self.tech_covs]
        self.execute_command(command)

    def execute_command(self, command):
        self.print_string("{}".format(' '.join(command)))
        if self.verbose:
            output = subprocess.check_call(command)
        else:
            output = subprocess.check_call(command,
                                           stdout=open(os.devnull, 'w'),
                                           stderr=subprocess.STDOUT)
        if output != 0:
            if not self.verbose:
                print("{}".format(' '.join(command)))
            print(output)
            exit()

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Technical covariates: {}".format(' '.join(self.tech_covs)))
        print("  > Interaction analyser: {}".format(self.eqtl_ia))
        print("  > Groups: {}".format(self.groups))
        print("  > Force: {}".format(self.force))
        print("  > Verbose: {}".format(self.verbose))
        print("")
