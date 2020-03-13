#!/usr/bin/env python3

"""
File:         eqtl_interaction_analyser_per_group.py
Created:      2020/03/06
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
import os

# Third party imports.

# Local application imports.


# Metadata.
__program__ = "eQTL Interaction Analyser (per group)"
__author__ = "M. Vochteloo"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class Main:
    """
    Main class of the program.
    """

    def __init__(self, indir, eqtl_ia):
        """
        Initializer method for the main class.

        :param indir: string, the folder with the input directories.
        :param eqtl_ia: string, the path to the interaction analyser.
        """
        self.indir = indir
        self.eqtl_ia = eqtl_ia

        self.outdir = os.path.join(os.getcwd(), 'output_w_corr')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        """
        Main method for the main class. Does all the work.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > eQTL Interaction Analyser: {}".format(self.eqtl_ia))
        print("  > Output directory: {}".format(self.outdir))
        print("")

        # Loop over the groups.
        for group_name in next(os.walk(self.indir))[1]:
            print("Working on: {}.".format(group_name))

            # create input / output directory.
            group_indir = os.path.join(self.indir, group_name)
            group_outdir = os.path.join(self.outdir, group_name)
            if not os.path.exists(group_outdir):
                os.makedirs(group_outdir)

            # create the eqtl filename.
            eqtl_file = os.path.join(group_indir, 'eqtl_table.txt.gz')

            # define the covariates to correct for.
            tech_covs = [
                "PCT_MRNA_BASES",
                "PCT_INTRONIC_BASES",
                "MEDIAN_3PRIME_BIAS",
                "PCT_USABLE_BASES",
                "PCT_INTERGENIC_BASES",
                "PCT_UTR_BASES",
                "PCT_READS_ALIGNED_IN_PAIRS",
                "PCT_CHIMERAS",
                "PF_READS_IMPROPER_PAIRS",
                "PF_HQ_ALIGNED_Q20_BASES",
                "PF_HQ_ALIGNED_BASES",
                "PCT_PF_READS_IMPROPER_PAIRS",
                "PF_READS_ALIGNED",
                "avg_mapped_read_length",
                "avg_input_read_length",
                "uniquely_mapped",
                "total_reads",
                "Total.Sequences_R1",
                "MDS1",
                "MDS2",
                "MDS3",
                "MDS4"
            ]

            # execute the program.
            command = 'java -jar {} -input {} -output {} -eqtls {} ' \
                      '-maxcov 1 ' \
                      '-noNormalization ' \
                      '-nnoCovNormalization ' \
                      '-cov {}'.format(self.eqtl_ia, group_indir,
                                       group_outdir, eqtl_file,
                                       ",".join(tech_covs))
            print("\t{}".format(command))
            os.system(command)


if __name__ == "__main__":
    INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen",
                         "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-03-03-interaction-analyser",
                         "step5-prepare-ia-inputs", "output_p_snp",
                         "groups")

    EQTL_IA = os.path.join(os.path.sep, "groups", "umcg-biogen",
                           "tmp03", "output",
                           "2019-11-06-FreezeTwoDotOne",
                           "2020-03-03-interaction-analyser",
                           "step6-interaction-analyser",
                           "eQTLInteractionAnalyser-1.2-SNAPSHOT-jar-with-dependencies.jar")

    # Start the program.
    MAIN = Main(indir=INDIR,
                eqtl_ia=EQTL_IA)

    MAIN.start()
