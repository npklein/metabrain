#!/usr/bin/env python3

"""
File:         decon_optimizer.py
Created:      2020/11/16
Last Changed: 2020/11/18
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

"""
./decon_optimizer.py -eq ../2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ge ../2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/create_matrices/genotype_table.txt.gz -al ../2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/create_matrices/genotype_alleles.txt.gz -ex ../2020-03-12-deconvolution/matrix_preparation/cis_new_output_log2/create_matrices/expression_table.txt.gz -cf data/deconvolution_table.txt.gz -de test_scripts/output/cis_cortex_deconvolutionResults_withFDR.txt.gz -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -sid rnaseq_id -cid cohort
"""

# Standard imports.

# Third party imports.

# Local application imports.
from src.main import Main
from src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Deconvolution Optimizer"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)
if __name__ == '__main__':
    # Get the command line arguments.
    CLA = CommandLineArguments(program=__program__,
                               version=__version__,
                               description=__description__)
    EQTL_PATH = CLA.get_argument('eqtl')
    GENOTYPE_PATH = CLA.get_argument('genotype')
    ALLELES_PATH = CLA.get_argument('alleles')
    EXPRESSION_PATH = CLA.get_argument('expression')
    CELL_FRACTIONS_PATH = CLA.get_argument('cell_fractions')
    DECON_PATH = CLA.get_argument('decon')
    SAMPLE_ANNOTATION_PATH = CLA.get_argument('sample_annotation')
    SAMPLE_ID = CLA.get_argument('sample_id')
    COHORT_ID = CLA.get_argument('cohort_id')
    ALPHA = CLA.get_argument('alpha')
    OUTDIR = CLA.get_argument('outdir')
    EXTENSIONS = CLA.get_argument('extension')

    # Start the program.
    PROGRAM = Main(eqtl_path=EQTL_PATH,
                   genotype_path=GENOTYPE_PATH,
                   alleles_path=ALLELES_PATH,
                   expression_path=EXPRESSION_PATH,
                   cell_fractions_path=CELL_FRACTIONS_PATH,
                   decon_path=DECON_PATH,
                   sample_annotation_path=SAMPLE_ANNOTATION_PATH,
                   sample_id=SAMPLE_ID,
                   cohort_id=COHORT_ID,
                   alpha=ALPHA,
                   outdir=OUTDIR,
                   extensions=EXTENSIONS)
    PROGRAM.start()