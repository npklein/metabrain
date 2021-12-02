#!/usr/bin/env python3

"""
File:         comnpare_cf_predictions.py
Created:      2021/11/22
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
__program__ = "Compare Cell Fraction Predictions"
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
./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/partial_deconvolution/PSYCHENCODE_PROFILE_CORTEX_EUR_TMM_LOG2/IHC_0CPM_LOG2/deconvolution_raw.txt.gz -transpose1 -n1 MetaBrain -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/Cell_fractions_Raw.txt.gz -n2 PsychENCODE_MyPredictions -o MetaBrain_vs_PsychENCODE_MyPredictions_CellFraction_comparison

./compare_cf_predictions.py -cf1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/partial_deconvolution/PSYCHENCODE_PROFILE_CORTEX_EUR_TPM_LOG2/NoENA_IHC_0CPM_LOG2/deconvolution_raw.txt.gz -transpose1 -n1 MetaBrain_TPM_NoENA -cf2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/PsychENCODE/Cell_fractions_Raw.txt.gz -n2 PsychENCODE_MyPredictions -o MetaBrain_TPM_NoENA_vs_PsychENCODE_MyPredictions_CellFraction_comparison

"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cell_fractions1_path = getattr(arguments, 'cell_fractions1')
        self.transpose1 = getattr(arguments, 'transpose1')
        self.name1 = getattr(arguments, 'name1')
        self.cell_fractions2_path = getattr(arguments, 'cell_fractions2')
        self.transpose2 = getattr(arguments, 'transpose2')
        self.name2 = getattr(arguments, 'name2')
        self.outfile = getattr(arguments, 'outfile')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

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
        parser.add_argument("-cf1",
                            "--cell_fractions1",
                            type=str,
                            required=True,
                            help="The path to the first cell fraction matrix")
        parser.add_argument("-transpose1",
                            action='store_true',
                            help="Transpose 1.")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first profile matrix")
        parser.add_argument("-cf2",
                            "--cell_fractions2",
                            type=str,
                            required=True,
                            help="The path to the second cell fraction matrix")
        parser.add_argument("-transpose2",
                            action='store_true',
                            help="Transpose 2.")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the second profile matrix")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=False,
                            default="deconvolution",
                            help="The name of the plot file.")

        return parser.parse_args()

    def start(self):
        print("Loading profile 1.")
        cf1_df = pd.read_csv(self.cell_fractions1_path, sep="\t", header=0, index_col=0)
        if self.transpose1:
            cf1_df = cf1_df.T
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cell_fractions1_path),
                                      cf1_df.shape))
        print(cf1_df)

        print("Loading profile 2.")
        cf2_df = pd.read_csv(self.cell_fractions2_path, sep="\t", header=0, index_col=0)
        if self.transpose2:
            cf2_df = cf2_df.T
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.cell_fractions2_path),
                                      cf2_df.shape))
        print(cf2_df)

        print("Subset overlap in samples")
        col_overlap = set(cf1_df.columns).intersection(set(cf2_df.columns))
        print("\tN-overlap: {}".format(len(col_overlap)))
        cf1_df = cf1_df.loc[:, col_overlap].copy()
        cf2_df = cf2_df.loc[:, col_overlap].copy()

        print("Correlating.")
        corr_df, pvalue_df = self.correlate(index_df=cf1_df.T, columns_df=cf2_df.T)
        signif_df = self.mask_non_significant(df=corr_df, pvalue_df=pvalue_df)
        self.plot_heatmap(df=signif_df,
                          annot_df=corr_df.round(2),
                          xlabel=self.name2,
                          ylabel=self.name1)

        print("Plotting.")
        cell_types = set(cf1_df.index).intersection(set(cf2_df.index))
        for ct in cell_types:
            plot_df = cf1_df.loc[[ct], :].T.merge(cf2_df.loc[[ct], :].T, left_index=True, right_index=True)
            plot_df.columns = ["x", "y"]

            self.plot_regplot(df=plot_df,
                              xlabel=self.name1,
                              ylabel=self.name2,
                              title=ct,
                              filename="{}_{}".format(self.outfile, ct))

    @staticmethod
    def correlate(index_df, columns_df, triangle=False):
        index_df_n_nan_values = index_df.shape[0] - index_df.isna().sum(axis=0)
        column_df_n_nan_values = columns_df.shape[0] - columns_df.isna().sum(axis=0)

        index_df_colnames = ["{} [N={:,}]".format(colname, index_df_n_nan_values[colname]) for colname in index_df.columns]
        column_df_colnames = ["{} [N={:,}]".format(colname, column_df_n_nan_values[colname]) for colname in columns_df.columns]

        corr_df = pd.DataFrame(np.nan, index=index_df_colnames, columns=column_df_colnames)
        pvalue_df = pd.DataFrame(np.nan, index=index_df_colnames, columns=column_df_colnames)

        for i, (index_column, index_colname) in enumerate(zip(index_df.columns, index_df_colnames)):
            for j, (column_column, column_colname) in enumerate(zip(columns_df.columns, column_df_colnames)):
                if triangle and i < j:
                    continue
                corr_data = pd.concat([index_df[index_column], columns_df[column_column]], axis=1)
                corr_data.dropna(inplace=True)

                coef = np.nan
                pvalue = np.nan
                if np.min(corr_data.std(axis=0)) > 0:
                    coef, pvalue = stats.pearsonr(corr_data.iloc[:, 1], corr_data.iloc[:, 0])

                corr_df.loc[index_colname, column_colname] = coef
                pvalue_df.loc[index_colname, column_colname] = pvalue

        return corr_df, pvalue_df

    @staticmethod
    def mask_non_significant(df, pvalue_df, a=0.05):
        signif_df = df.copy()
        for i in range(signif_df.shape[0]):
            for j in range(signif_df.shape[1]):
                if np.isnan(pvalue_df.iloc[i, j]) or pvalue_df.iloc[i, j] >= a:
                    signif_df.iloc[i, j] = 0

        return signif_df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="", vmin=-1, vmax=1):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * df.shape[1] + 10, 1 * df.shape[0] + 10),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        annot_df.fillna("", inplace=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:
                sns.heatmap(df, cmap=cmap, vmin=vmin, vmax=vmax, center=0,
                            square=True, annot=annot_df, fmt='',
                            cbar=False, annot_kws={"size": 14, "color": "#000000"},
                            ax=ax)

                plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
                plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

                ax.set_xlabel(xlabel, fontsize=14)
                ax.xaxis.set_label_position('top')

                ax.set_ylabel(ylabel, fontsize=14)
                ax.yaxis.set_label_position('right')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > 1:
                col_index = 0
                row_index += 1

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_corr_heatmap_Pearson.png".format(self.outfile)))
        plt.close()

    def plot_regplot(self, df, x="x", y="y", xlabel="", ylabel="", title="",
                     filename="plot"):

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # Plot.
        sns.regplot(x=x, y=y, data=df, ci=95,
                    scatter_kws={'facecolors': "#000000",
                                 'linewidth': 0,
                                 'alpha': 0.75},
                    line_kws={"color": "#b22222"},
                    ax=ax
                    )

        # Regression.
        coef_str = "NA"
        if np.std(df[y]) != 0 and np.std(df[x]) != 0:
            coef, _ = stats.pearsonr(df[y], df[x])
            coef_str = "{:.2f}".format(coef)

        # Add the text.
        ax.annotate(
            'N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=0.75,
            fontsize=14,
            fontweight='bold')
        ax.annotate(
            'r = {}'.format(coef_str),
            xy=(0.03, 0.9),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=0.75,
            fontsize=14,
            fontweight='bold')

        ax.set_title(title,
                     fontsize=18,
                     fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}_regplot.png".format(filename)))
        print("\tSaved '{}'".format(filename))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
