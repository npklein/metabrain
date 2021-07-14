#!/usr/bin/env python3

"""
File:         compare_deconvolution_matrices.py
Created:      2020/09/03
Last Changed: 2021/07/14
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
./compare_deconvolution_matrices.py -d1 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n1 Niek -d2 ../2021-06-24-decon-QTL/CortexEUR-cis/deconvolutionResults.csv -n2 Martijn -o Niek_vs_Martijn

./compare_deconvolution_matrices.py -d1 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n1 old_pre_processing -d2 ../2021-06-24-decon-QTL/CortexEUR-cis-OldCellCounts/deconvolutionResults.csv -n2 new_pre_processing -o OldPreProcessing_vs_NewPreProcessing

./compare_deconvolution_matrices.py -d1 ../2021-06-24-decon-QTL/CortexEUR-cis-OldCellCounts/deconvolutionResults.csv -n1 old_cell-counts -d2 ../2021-06-24-decon-QTL/CortexEUR-cis/deconvolutionResults.csv -n2 new_cell_counts -o OldCellCounts_vs_NewCellCounts

./compare_deconvolution_matrices.py -d1 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n1 Decon-eQTL -d2 ../decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/output/deconvolutionResults.txt.gz -n2 Martijn -o DeconeQTL_vs_Martijn

./compare_deconvolution_matrices.py -d1 ../decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-woPermFDR-Old/deconvolutionResults.txt.gz -n1 OLD -d2 ../decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-woPermFDR-New/deconvolutionResults.txt.gz -n2 NEW -o PythonCode_Old_vs_New

./compare_deconvolution_matrices.py -d1 ../decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-woPermFDR-New-OldCC/deconvolutionResults.txt.gz -n1 OLD_cc -d2 ../decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-woPermFDR-New/deconvolutionResults.txt.gz -n2 NEW_cc -o PythonCode_OldCC_vs_NewCC
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon1_path = getattr(arguments, 'decon1')
        self.decon1_name = getattr(arguments, 'name1')
        self.decon2_path = getattr(arguments, 'decon2')
        self.decon2_name = getattr(arguments, 'name2')
        self.log10_transform = getattr(arguments, 'log10')
        self.outfile = getattr(arguments, 'outfile')
        self.nrows = None

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.accent_colormap = {
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Average": "#E8E8E8"
        }
        self.point_colormap = {
            "no signif": "#808080",
            "x signif": "#0072B2",
            "y signif": "#D55E00",
            "both signif": "#009E73"
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
        parser.add_argument("-log10",
                            action='store_true',
                            help="Log10 transform the profile values."
                                 " Default: False.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=False,
                            default="deconvolution",
                            help="The name of the plot file.")

        return parser.parse_args()

    def start(self):
        print("Loading deconvolution matrix 1.")
        decon1_df = pd.read_csv(self.decon1_path, sep="\t", header=0,
                                index_col=0, nrows=self.nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon1_path),
                                      decon1_df.shape))
        decon1_pval_df, decon1_fdr_df = self.bh_correct(decon1_df)
        # decon2_pval_df, decon2_fdr_df = self.bh_correct2(decon1_df)

        print("Loading deconvolution matrix 2.")
        decon2_df = pd.read_csv(self.decon2_path, sep="\t", header=0,
                                index_col=0, nrows=self.nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon2_path),
                                      decon2_df.shape))
        decon2_pval_df, decon2_fdr_df = self.bh_correct(decon2_df)

        for decon1_df, decon2_df, appendix in ([decon1_pval_df, decon2_pval_df, "pvalue"], [decon1_fdr_df, decon2_fdr_df, "FDR"]):
            print("Plotting {}".format(appendix))

            # Change columns.
            decon1_df.columns = [x.split("_")[1] for x in decon1_df.columns]
            decon2_df.columns = [x.split("_")[1] for x in decon2_df.columns]

            # Filter on overlapping rows and columns.
            print("\tFiltering matrices.")
            row_overlap = set(decon1_df.index).intersection(set(decon2_df.index))
            col_overlap = set(decon1_df.columns).intersection(set(decon2_df.columns))

            decon1_df = decon1_df.loc[row_overlap, col_overlap]
            decon2_df = decon2_df.loc[row_overlap, col_overlap]

            print("\t  Deconvolution matrix 1: {}".format(decon1_df.shape))
            print("\t  Deconvolution matrix 2: {}".format(decon2_df.shape))

            if decon1_df.shape != decon2_df.shape:
                print("\t  Shape's are not identical.")
                exit()

            if decon1_df.shape[0] > decon1_df.shape[1]:
                decon1_df = decon1_df.T
                decon2_df = decon2_df.T

            print("\tPlotting.")
            #self.plot_violin_comparison(decon1_df, decon2_df, self.outfile + "_" + appendix)
            self.plot_regression_comparison(decon1_df, decon2_df, self.outfile + "_" + appendix)

    @staticmethod
    def bh_correct(input_df):
        pvalue_cols = []
        for col in input_df.columns:
            if col.endswith("_pvalue"):
                pvalue_cols.append(col)
        pvalue_df = input_df.loc[:, pvalue_cols]

        fdr_data = []
        indices = []
        for col in pvalue_cols:
            fdr_data.append(multitest.multipletests(input_df.loc[:, col], method='fdr_bh')[1])
            indices.append(col.replace("_pvalue", ""))
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=input_df.index)

        return pvalue_df, fdr_df.T

    @staticmethod
    def bh_correct2(input_df):
        pvalue_cols = []
        for col in input_df.columns:
            if col.endswith("_pvalue"):
                pvalue_cols.append(col)
        pvalue_df = input_df.loc[:, pvalue_cols]

        pvalue_df_m = pvalue_df.copy()
        pvalue_df_m.reset_index(drop=False, inplace=True)
        pvalue_df_m = pvalue_df_m.melt(id_vars=["index"])
        pvalue_df_m["FDR"] = multitest.multipletests(pvalue_df_m.loc[:, "value"], method='fdr_bh')[1]
        fdr_df = pvalue_df_m.pivot(index='index', columns='variable', values='value')
        fdr_df.columns = [x.replace("_pvalue", "") for x in fdr_df.columns]

        return pvalue_df, fdr_df

    def plot_violin_comparison(self, decon1_df, decon2_df, filename):
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
        fig.savefig(os.path.join(self.outdir, "{}_violin_comparison.png".format(filename)))
        plt.close()

    def plot_regression_comparison(self, decon1_df, decon2_df, filename):
        # Calculate number of rows / columns.
        celltypes = list(decon1_df.index)
        celltypes.sort()
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

            # Add color
            df["hue"] = self.point_colormap["no signif"]
            df.loc[(df[x_label] < 0.05) & (df[y_label] >= 0.05), "hue"] = self.point_colormap["x signif"]
            df.loc[(df[x_label] >= 0.05) & (df[y_label] < 0.05), "hue"] = self.point_colormap["y signif"]
            df.loc[(df[x_label] < 0.05) & (df[y_label] < 0.05), "hue"] = self.point_colormap["both signif"]
            counts = df["hue"].value_counts()
            for value in self.point_colormap.values():
                if value not in counts.index:
                    counts[value] = 0

            coef, _ = stats.spearmanr(df[x_label], df[y_label])
            coef_str = "{:.2f}".format(coef)

            accent_color = self.accent_colormap[celltype]

            # Convert to log10 scale.
            if self.log10_transform:
                df.loc[df[x_label] == 0, x_label] = 2.2250738585072014e-308
                df.loc[df[y_label] == 0, y_label] = 2.2250738585072014e-308
                df[[x_label, y_label]] = np.log10(df[[x_label, y_label]])

            # Creat the subplot.
            ax = fig.add_subplot(grid[row_index, col_index])
            sns.despine(fig=fig, ax=ax)

            # Plot.
            g = sns.regplot(x=x_label,
                            y=y_label,
                            data=df,
                            scatter_kws={'facecolors': df["hue"],
                                         'linewidth': 0,
                                         'alpha': 0.5},
                            line_kws={"color": accent_color},
                            ax=ax
                            )

            # Add the text.
            ax.annotate(
                'r = {}'.format(coef_str),
                xy=(0.03, 0.94),
                xycoords=ax.transAxes,
                color="#404040",
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'total N = {}'.format(df.shape[0]),
                xy=(0.03, 0.9),
                xycoords=ax.transAxes,
                color="#404040",
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {}'.format(counts[self.point_colormap["both signif"]]),
                xy=(0.03, 0.86),
                xycoords=ax.transAxes,
                color=self.point_colormap["both signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {}'.format(counts[self.point_colormap["x signif"]]),
                xy=(0.03, 0.82),
                xycoords=ax.transAxes,
                color=self.point_colormap["x signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {}'.format(counts[self.point_colormap["y signif"]]),
                xy=(0.03, 0.78),
                xycoords=ax.transAxes,
                color=self.point_colormap["y signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {}'.format(counts[self.point_colormap["no signif"]]),
                xy=(0.03, 0.74),
                xycoords=ax.transAxes,
                color=self.point_colormap["no signif"],
                fontsize=18,
                fontweight='bold')

            ax.text(0.5, 1.05,
                    celltype,
                    fontsize=26, weight='bold', ha='center', va='bottom',
                    color=accent_color,
                    transform=ax.transAxes)

            ax.set_xlabel(" ".join(x_label.split("_")),
                          fontsize=18,
                          fontweight='bold')

            ax.set_ylabel(" ".join(y_label.split("_")),
                          fontsize=18,
                          fontweight='bold')

            signif_line = 0.05
            if self.log10_transform:
                signif_line = np.log10(0.05)
            ax.axhline(signif_line, ls='--', color="#000000", zorder=-1)
            ax.axvline(signif_line, ls='--', color="#000000", zorder=-1)

            min_value = np.min((df[[x_label, y_label]].min(axis=0).min() * 1.1, -0.1))
            max_value = np.max((df[[x_label, y_label]].max(axis=0).max() * 1.1, 1.1))

            ax.set_xlim(min_value, max_value)
            ax.set_ylim(min_value, max_value)
            ax.plot([min_value, max_value], [min_value, max_value], ls="--", c="#000000")

            # Increment indices.
            col_index += 1
            if col_index == ncols:
                col_index = 0
                row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_regression_comparison.png".format(filename)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
