#!/usr/bin/env python3

"""
File:         compare_deconvolution_matrices.py
Created:      2020/09/03
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
import math
import os

# Third party imports.
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Compare Deconvolution Matrices"
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
./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_perIeQTLFDR.txt.gz -n1 per_ieqtl -d2 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_CombinedPermFDR.txt.gz -n2 combined

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_perIeQTLFDR.txt.gz -n1 per_ieqtl -d2 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_CombinedPermFDR.txt.gz -n2 combined

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation/decon_output/real/deconvolutionResults.csv -n1 p_value -d2 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_CombinedPermFDR.txt.gz -n2 combined_FDR

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation/decon_output/real/deconvolutionResults.csv -n1 p_value -d2 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_perIeQTLFDR.txt.gz -n2 single_FDR

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation/decon_output/real/deconvolutionResults.csv -n1 p_value -d2 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_perIeQTLFDR.txt.gz -n2 single_FDR

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation/decon_output/real/deconvolutionResults.csv -n1 p_value -d2 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n2 BH_FDR

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_perIeQTLFDR.txt.gz -n1 single_FDR -d2 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n2 BH_FDR

./compare_deconvolution_matrices.py -d1 ../paper_scripts/decon_eqtl_permutation_per_cohort/deconvolutionResults_1000perm_CombinedPermFDR.txt.gz -n1 combined_FDR -d2 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n2 BH_FDR

./compare_deconvolution_matrices.py -d1 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n1 Niek -d2 ../2021-06-24-decon-QTL/cortex_eur_cis/decon_cis_cortex_eur_out/deconvolutionResults.csv -n2 Martijn
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon1_path = getattr(arguments, 'decon1')
        self.decon1_name = getattr(arguments, 'name1')
        self.decon2_path = getattr(arguments, 'decon2')
        self.decon2_name = getattr(arguments, 'name2')
        self.nrows = None

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.colormap = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Average": "#E8E8E8"
        }

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-d1",
                            "--decon1",
                            type=str,
                            required=True,
                            help="The path to the first deconvolution matrix")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first deconvolution matrix")
        parser.add_argument("-d2",
                            "--decon2",
                            type=str,
                            required=True,
                            help="The path to the second deconvolution matrix")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the first deconvolution matrix")

        return parser.parse_args()

    def start(self):
        print("Loading deconvolution matrix 1.")
        decon1_df = pd.read_csv(self.decon1_path, sep="\t", header=0,
                                index_col=0, nrows=self.nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon1_path),
                                      decon1_df.shape))
        decon1_df, _ = self.bh_correct(decon1_df)

        print("Loading deconvolution matrix 2.")
        decon2_df = pd.read_csv(self.decon2_path, sep="\t", header=0,
                                index_col=0, nrows=self.nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon2_path),
                                      decon2_df.shape))
        decon2_df, _ = self.bh_correct(decon2_df)

        print(decon1_df)
        print(decon2_df)

        # Change columns.
        decon1_df.columns = [x.split("_")[1] for x in decon1_df.columns]
        decon2_df.columns = [x.split("_")[1] for x in decon2_df.columns]

        # Filter on overlapping rows and columns.
        print("Filtering matrices.")
        row_overlap = set(decon1_df.index).intersection(set(decon2_df.index))
        col_overlap = set(decon1_df.columns).intersection(set(decon2_df.columns))

        decon1_df = decon1_df.loc[row_overlap, col_overlap]
        decon2_df = decon2_df.loc[row_overlap, col_overlap]

        print("\tDeconvolution matrix 1: {}".format(decon1_df.shape))
        print("\tDeconvolution matrix 2: {}".format(decon2_df.shape))

        if decon1_df.shape != decon2_df.shape:
            print("Shape's are not identical.")
            exit()

        if decon1_df.shape[0] > decon1_df.shape[1]:
            decon1_df = decon1_df.T
            decon2_df = decon2_df.T

        #self.plot_violin_comparison(decon1_df, decon2_df)
        self.plot_regression_comparison(decon1_df, decon2_df)

    @staticmethod
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        pval_data = []
        fdr_data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                pval_data.append(df.loc[:, col])
                fdr_data.append(multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        pval_df = pd.DataFrame(pval_data, index=indices, columns=df.index)
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=df.index)

        return pval_df.T, fdr_df.T

    def plot_violin_comparison(self, decon1_df, decon2_df):
        df1 = decon1_df.copy()
        df2 = decon2_df.copy()

        df1 = df1.T.melt()
        df2 = df2.T.melt()

        df1["name"] = self.decon1_name
        df2["name"] = self.decon2_name

        df = pd.concat([df1, df2], axis=0)
        del df1, df2

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.violinplot(x="variable", y="value", hue="name",
                       data=df, ax=ax)
        ax.set_xlabel("",
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("cell count %",
                      fontsize=14,
                      fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "deconvolution_violin_comparison.png"))
        plt.close()

    def plot_regression_comparison(self, decon1_df, decon2_df):
        # Plot.
        print("Plotting")

        # Calculate number of rows / columns.
        celltypes = list(decon1_df.index)
        ncols = 3
        nrows = math.ceil(len(celltypes) / 3)

        sns.set(rc={'figure.figsize': (12 * ncols, 9 * nrows)})
        sns.set_style("ticks")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=ncols,
                                nrows=nrows)

        row_index = 0
        col_index = 0

        for celltype in celltypes:
            x_label = "{}_{}".format(self.decon1_name, celltype)
            y_label = "{}_{}".format(self.decon2_name, celltype)
            df = pd.DataFrame({x_label: decon1_df.loc[celltype, :],
                               y_label: decon2_df.loc[celltype, :]})

            # remove not signif.
            #df = df.loc[df.min(axis=1) < 0.05, :]
            df = np.log10(df + 1)

            coef, p_value = stats.spearmanr(df[x_label], df[y_label])
            coef_str = "{:.2f}".format(coef)
            p_str = "p = {:.2e}".format(p_value)

            color = self.colormap[celltype]

            # Creat the subplot.
            ax = fig.add_subplot(grid[row_index, col_index])
            sns.despine(fig=fig, ax=ax)

            # Plot.
            g = sns.regplot(x=x_label,
                            y=y_label,
                            data=df,
                            scatter_kws={'facecolors': '#000000',
                                         'linewidth': 0,
                                         'alpha': 0.5},
                            line_kws={"color": color},
                            ax=ax
                            )

            # Add the text.
            ax.annotate(
                'r = {} [{}]'.format(coef_str, p_str),
                xy=(0.03, 0.94),
                xycoords=ax.transAxes,
                color=color,
                fontsize=18,
                fontweight='bold')

            ax.text(0.5, 1.05,
                    celltype,
                    fontsize=26, weight='bold', ha='center', va='bottom',
                    color=color,
                    transform=ax.transAxes)

            ax.set_xlabel(" ".join(x_label.split("_")),
                          fontsize=18,
                          fontweight='bold')

            ax.set_ylabel(" ".join(y_label.split("_")),
                          fontsize=18,
                          fontweight='bold')

            # ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
            # ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
            #
            # ax.set_xlim(-0.1, 1.1)
            # ax.set_ylim(-0.1, 1.1)
            # ax.plot([-0.1, 1.1], [-0.1, 1.1], ls="--", c=".3")

            # Increment indices.
            col_index += 1
            if col_index == ncols:
                col_index = 0
                row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "deconvolution_regression_comparison.png"))
        plt.close()

if __name__ == '__main__':
    m = main()
    m.start()
