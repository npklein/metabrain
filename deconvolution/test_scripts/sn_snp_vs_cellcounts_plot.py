#!/usr/bin/env python3

"""
File:         sn_snp_vs_cellcounts_plot.py
Created:      2020/11/27
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
import argparse
import gzip
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "SN SNP vs CellCounts Plot"
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

"""
Syntax:
./sn_snp_vs_cellcounts_plot.py -eq /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/trans_100Perm/EX/eQTLProbesFDR0.05-ProbeLevel.txt.gz -ge /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_trans/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/EX_expression.txt -cc /groups/umcg-biogen/tmp03/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/cell_counts.txt -gec /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/ROSMAP-scRNAseq-genometoexpressioncoupling-withBulkID.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.sn_expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.gec_path = getattr(arguments, 'geno_expr_coupling')
        self.extensions = getattr(arguments, 'extension')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), "sn_snp_vs_cellcounts_plot")
        self.palette = {
            "In": "#56B4E9",
            "Ex": "#0072B2",
            "Oli": "#009E73",
            "End": "#CC79A7",
            "Mic": "#E69F00",
            "Ast": "#D55E00",
            "Per": "#808080",
            "Opc": "#F0E442"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix")
        parser.add_argument("-gec",
                            "--geno_expr_coupling",
                            type=str,
                            required=False,
                            default=None,
                            help="A couple file for genotype - expression")
        parser.add_argument("-e",
                            "--extension",
                            type=str,
                            choices=["png", "pdf"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        eqtl_df = self.load_file(self.eqtl_path, index_col=None)
        print(eqtl_df)
        geno_df = self.load_file(self.geno_path)
        print(geno_df)
        expr_df = self.load_file(self.expr_path)
        print(expr_df)
        cc_df = self.load_file(self.cc_path)
        print(cc_df)
        gec_df = None
        if self.gec_path is not None:
            gec_df = self.load_file(self.gec_path, index_col=None)
            print(gec_df)

        print("### Step2 ###")
        print("Filter data")

        if gec_df is not None:
            geno_df = geno_df.loc[:, gec_df.loc[:, "bulk-expr"]].copy()
            geno_df.columns = [str(x) for x in gec_df.loc[:, "sn-expr"].values]
            geno_df.drop_duplicates(inplace=True)
            geno_df.index.name = "-"
            print(geno_df)

        print("### Step3 ###")
        print("Visualise snp vs cc")
        plotted_snps = set()
        for i, (index, row) in enumerate(eqtl_df.iterrows()):
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]

            if snp_name in plotted_snps:
                continue

            plotted_snps.add(snp_name)

            print("\tWorking on: {}\t{}\t{}".format(snp_name, probe_name, hgnc_name))

            # Get the genotype / expression data.
            genotype = geno_df.loc[[snp_name], :].copy()
            genotype.index = ["genotype"]
            expression = expr_df.loc[[probe_name], :].copy()
            expression.index = ["expression"]
            data = genotype.T.merge(expression.T, left_index=True, right_index=True).merge(cc_df.T, left_index=True, right_index=True)

            # Remove missing values.
            data = data.loc[(data['genotype'] >= 0.0) &
                            (data['genotype'] <= 2.0), :]

            for ct in cc_df.index:
                if ct != "Ex":
                    continue

                cell_data = data[["genotype", ct]]
                cell_data.columns = ["x", "y"]

                self.plot(snp_name, ct, cell_data)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def plot(self, snp_name, ct, df):
        # Calculate the correlation.
        coef, p = stats.spearmanr(df["x"], df["y"])

        # Prepare the figure.
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot the scatter / box plot.
        sns.regplot(x="x", y="y", data=df,
                    scatter=False,
                    line_kws={"color": self.palette[ct]},
                    ax=ax
                    )
        sns.boxplot(x="x", y="y", data=df,
                    palette={0.0: "#D55E00", 1.0: "#0072B2", 2.0: "#E69F00"},
                    showfliers=False,
                    zorder=1,
                    # boxprops=dict(alpha=.3),
                    ax=ax)

        # Set the other aesthetics.
        ax.text(0.5, 1.06,
                "{} vs {} cell count".format(snp_name, ct),
                fontsize=22, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'r = {:.2f} [p = {:.2e}]'.format(coef, p),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_ylabel('{} cell count'.format(ct),
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel('SNP {}'.format(snp_name),
                      fontsize=14,
                      fontweight='bold')

        # Safe the plot.
        for extension in self.extensions:
            filename = "{}_{}.{}".format(snp_name, ct, extension)
            print("\t\tSaving plot: {}".format(filename))
            fig.savefig(os.path.join(self.outdir, filename))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Gentoype expression coupling path: {}".format(self.gec_path))
        print("  > Extension: {}".format(self.extensions))
        print("  > Output directory: {}".format(self.outdir))
        print("")



if __name__ == '__main__':
    m = main()
    m.start()
