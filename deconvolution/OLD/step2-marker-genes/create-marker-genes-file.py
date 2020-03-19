#!/usr/bin/env python3

"""
File:         create_marker_genes_file.py
Created:      2020/03/11
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2019 M.Vochteloo

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
import numpy as np
import pandas as pd

# Local application imports.


# Metadata.
__program__ = "Create Marker Genes File"
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

    def __init__(self, eqtl_file, marker_dict, expr_file, snp_translate):
        """
        The initializer for the main class.

        :param eqtl_file: string, the file containing the eQTL effects of
                          interest + also acts as translation file.
        :param marker_dict: dict, list of cell types: marker genes.
        :param expr_file: string, the expression data input file.
        :param snp_translate: string, SNP translate file.
        """
        self.eqtl_file = eqtl_file
        self.marker_dict = marker_dict
        self.expr_file = expr_file
        self.snp_translate = snp_translate
        self.outdir = os.path.join(os.getcwd(), 'output')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Print arguments.
        print("Arguments:")
        print("  > eQTL input file: {}".format(self.eqtl_file))
        print("  > Marker dict: {}".format(self.marker_dict))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > SNP translate file: {}".format(self.snp_translate))
        print("  > Output directory: {}".format(self.outdir))
        print("")

        # Load the eQTL file.
        print("Loading eQTL matrix.")
        eqtl_df = pd.read_csv(self.eqtl_file, sep="\t", header=0)
        print("\tShape: {}".format(eqtl_df.shape))

        # Load the expression matrix file.
        print("Loading expression matrix.")
        expr_df = pd.read_csv(self.expr_file, sep="\t", header=0, nrows=nrows)
        print("\tShape: {}".format(expr_df.shape))

        # Load the SNP translate file.
        print("Loading SNP translate matrix.")
        translate_df = pd.read_csv(self.snp_translate, sep="\t", header=0)
        print("\tShape: {}".format(translate_df.shape))

        # Add the masked name to the eQTL df.
        eqtl_df = pd.merge(eqtl_df, translate_df, on=['SNPName', 'ProbeName'])
        eqtl_df.index = eqtl_df["HGNCName"]
        eqtl_df.index.name = "HGNCName"

        # Loop over the celltypes / marker genes.
        eqtl_outlist = []
        cov_out = pd.DataFrame(columns=expr_df.columns.to_list())
        for celltype, marker_genes in self.marker_dict.items():
            for marker_gene in marker_genes:
                if marker_gene in eqtl_df.index:
                    # Get the info.
                    info = eqtl_df.loc[[marker_gene], ["SNPName", "ProbeName", "masked"]]
                    for i, (index, (snp_name, probe_name, masked_name)) in enumerate(info.iterrows()):
                        eqtl_outlist.append([celltype, marker_gene, snp_name, probe_name, masked_name])

                        if i == 0:
                            # Create a new row for the dataframe.
                            # probe_name = "ENSG00000000003.15"
                            expression = expr_df.loc[expr_df["-"] == probe_name, :].copy()
                            expression.drop(["-"], axis=1, inplace=True)
                            expression.index = ["{}_{}".format(celltype, marker_gene)]
                            cov_out = pd.concat([expression, cov_out])

        # write marker gene covariate outfile.
        cov_out.index.name = "-"
        cov_out.to_csv(os.path.join(self.outdir, "marker_genes.txt.gz"),
                       sep="\t", compression='gzip')

        # write marker gene eQTLs outfile.
        eqtl_out = pd.DataFrame(eqtl_outlist, columns=["CellType", "HGNCName",
                                        "SNPName", "ProbeName", "MaskedName"])
        eqtl_out.index.name = "-"
        eqtl_out.to_csv(os.path.join(self.outdir, "marker_genes_eQTLs.txt.gz"),
                        sep="\t", compression='gzip')


if __name__ == "__main__":
    # Step 1 data.
    EQTLS = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-03-03-interaction-analyser",
                         "step1-combine-eQTLprobe-files", "output",
                         "eQTLProbesFDR0.05-ProbeLevel_combined.txt.gz")

    # Other data sources.
    EXPRESSION = os.path.join(os.path.sep, "groups", "umcg-biogen",
                              "tmp03", "output",
                              "2019-11-06-FreezeTwoDotOne",
                              "2020-01-31-expression-tables",
                              "2020-02-05-step6-covariate-removal",
                              "2020-02-17-step5-remove-covariates-per-dataset",
                              "output-cortex",
                              "MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemovedOLS.txt.gz"
                              )

    SNP_TRANSLATE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-03-03-interaction-analyser",
                         "step4-prepare-matrices", "output", "masked",
                         "SNP_translate_table.txt.gz")

    # Define marker genes.
    MARKERS = {"neurons": ["NPY", "VIP", "CCK", "CALB2", "KIT", "DCX", "SNAP25",
                           "UCHL1", "INA", "GAP43", "SST", "ERBB4", "SYT10",
                           "IGF1", "ROBO2", "RELN", "STMN2", "ELAVL3"],
               "oligodendrocytes": ["MOG", "MBP", "MAG", "CNP", "CA2", "PLP1",
                                    "OLIG1", "MOBP", "OMG", "ASPA", "SOX10",
                                    "GPR37", "OLIG2"],
               "endothelialcells": ["CD34", "VWF", "ENG", "KDR", "CXCL12",
                                    "FLT1", "SHE", "CAV1"],
               "microglia": ["AIF1", "CCL2", "CX3CR1", "TLR2", "CD14", "CD68",
                             "IL1A"],
               "astrocytes": ["GFAP", "AQP4", "EGFR", "SLC1A3", "SLC1A2",
                              "PAX6", "FABP7", "GJA1", "GJB6", "CLU", "MLC1"]}

    # Start the program.
    MAIN = Main(eqtl_file=EQTLS,
                marker_dict=MARKERS,
                expr_file=EXPRESSION,
                snp_translate=SNP_TRANSLATE)
    MAIN.start()
