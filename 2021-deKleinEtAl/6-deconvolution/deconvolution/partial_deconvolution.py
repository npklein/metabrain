#!/usr/bin/env python3

"""
File:         partial_deconvolution.py
Created:      2020/06/29
Last Changed: 2020/12/15
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

# Third party imports.

# Local application imports.
from partial_deconvolution.src.main import Main
from partial_deconvolution.src.settings import Settings
from partial_deconvolution.src.cmd_line_arguments import CommandLineArguments

# Metadata
__program__ = "Partial Deconvolution"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax:
./partial_deconvolution.py -d /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt -t /groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/IHC_counts.txt.gz -o ALL_TMM_LOG2 -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -sid rnaseq_id -cid cohort -os IHC_0CPM_LOG2_FILTERED_CC_TEST -cohort_corr -log2 -sum_to_one -visualise -e pdf

./partial_deconvolution.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-02-04-step5-center-scale/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/CellMap_brain_CNS7_avgCPM.txt -t /groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/AMP-AD/IHC_counts.txt.gz -o NEW_PROFILE_ALL_TMM_LOG2 -sa /groups/umcg-biogen/tmp01/output/2020-11-10-DeconOptimizer/2020-11-10-decon-optimizer/data/2021-04-28.brain.phenotypes.CortexEUR.txt.gz -sid rnaseq_id -cid MetaBrain_cohort -os IHC_0CPM_LOG2_FILTERED_CC -cohort_corr -log2 -sum_to_one -visualise -e pdf
"""


if __name__ == '__main__':
    # Get the command line arguments.
    CLA = CommandLineArguments(program=__program__,
                               version=__version__,
                               description=__description__)
    SETTINGS = Settings(data_path=CLA.get_argument("data"),
                        signature_path=CLA.get_argument("signature"),
                        translate_path=CLA.get_argument("translate"),
                        ground_truth_path=CLA.get_argument("ground_truth"),
                        sample_annotation_path=CLA.get_argument("sample_annotation"),
                        sample_id=CLA.get_argument("sample_id"),
                        sample_filter_path=CLA.get_argument("sample_filter"),
                        cohort_id=CLA.get_argument("cohort_id"),
                        cohort_filter=CLA.get_argument("cohort_filter"),
                        annotation_id=CLA.get_argument("annotation_id"),
                        annotation_filter=CLA.get_argument("annotation_filter"),
                        min_expr=CLA.get_argument("min_expr"),
                        cohort_corr=CLA.get_argument("cohort_corr"),
                        normalize=CLA.get_argument("normalize"),
                        zscore=CLA.get_argument("zscore"),
                        log2=CLA.get_argument("log2"),
                        decon_method=CLA.get_argument("decon_method"),
                        sum_to_one=CLA.get_argument("sum_to_one"),
                        extension=CLA.get_argument("extension")
                        )

    # Start the program.
    PROGRAM = Main(settings=SETTINGS,
                   outdir=CLA.get_argument("outdir"),
                   outsubdir=CLA.get_argument("outsubdir"),
                   visualise=CLA.get_argument("visualise"),
                   plot_ids=CLA.get_argument("plot_id"),)
    PROGRAM.start()
