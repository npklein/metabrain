#!/usr/bin/env python3

"""
File:         compare_zscores.py
Created:      2020/11/04
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
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats


# Local application imports.

# Metadata
__program__ = "Compare Zscores"
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


class main():
    def __init__(self):
        self.bulk_infolder = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-EUR/Iteration1/"
        self.bulk_filename = "eQTLProbesFDR0.05-ProbeLevel.txt.gz"
        self.decon_infile = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/decon_out/deconvolutionResults_withFDR.txt.gz"
        self.sn_infolder = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/cis/"
        self.sn_filename = "eQTLsFDR-ProbeLevel.txt.gz"
        self.cell_types = [("AST", "Astrocyte", "#D55E00"),
                           ("END", "EndothelialCell", "#CC79A7"),
                           ("EX", "Neuron", "#0072B2"),
                           ("IN", "Neuron", "#0072B2"),
                           ("MIC", "Macrophage", "#E69F00"),
                           ("OLI", "Oligodendrocyte", "#009E73")]
        self.extension = "png"

        self.outdir = os.path.join(str(Path(__file__).parent.parent),
                                   'plots')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        print("Loading data")

        bulk_df = pd.read_csv(os.path.join(self.bulk_infolder, self.bulk_filename),
                              sep="\t",
                              header=0,
                              index_col=None)
        print("\tBulk eQTL data frame: {}".format(bulk_df.shape))
        print(bulk_df)

        decon_df = pd.read_csv(self.decon_infile,
                               sep="\t",
                               header=0,
                               index_col=0,
                               nrows=bulk_df.shape[0])
        print("\tDeconvolution results data frame: {}".format(decon_df.shape))
        print(decon_df)
        print(list(decon_df.index))
        exit()

        print("Preprocessing deconvolution data frame")
        probe_names = []
        snp_names = []
        for index in decon_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        decon_df["ProbeName"] = probe_names
        decon_df["SNPName"] = snp_names
        decon_df.reset_index(drop=True, inplace=True)
        print(decon_df)

        print("Merging results")
        bulk_decon_df = bulk_df.merge(decon_df,
                                      left_on=["SNPName", "ProbeName"],
                                      right_on=["SNPName", "ProbeName"])
        del bulk_df
        bulk_decon_df.index = bulk_decon_df["SNPName"] + "_" + bulk_decon_df["ProbeName"]
        print(bulk_decon_df)

        print("Visualizing")

        sns.set(rc={'figure.figsize': (len(self.cell_types)*8, 18)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=3, ncols=len(self.cell_types))

        for col_index, (sn_ct, b_ct, color) in enumerate(self.cell_types):

            print("\tLoading single-nucleus results")
            sn_df = pd.read_csv(os.path.join(self.sn_infolder, sn_ct, self.sn_filename),
                                sep="\t",
                                header=0,
                                index_col=0)
            sn_df.index = sn_df["SNPName"] + "_" + sn_df["ProbeName"]
            print("\t\tSingle-nucleus eQTL data frame: {}".format(sn_df.shape))

            overlap = set(sn_df.index).intersection(set(bulk_decon_df.index))
            print("\t{} overlapping keys".format(len(overlap)))

            sn_df = sn_df.loc[overlap, :].copy()
            print("\t\tSingle-nucleus eQTL data frame: {}".format(sn_df.shape))
            ct_bulk_decon_df = bulk_decon_df.loc[overlap, :].copy()
            print("\t\t{} bulk deconvolution data frame: {}".format(b_ct, ct_bulk_decon_df.shape))

            print("\tComparing alleles")
            flip_mask = []
            for i in range(len(overlap)):
                if sn_df.iloc[i, :]["AlleleAssessed"] != \
                        ct_bulk_decon_df.iloc[i, :]["AlleleAssessed"]:
                    flip_mask.append(-1)
                else:
                    flip_mask.append(1)

            sn_df["flip"] = flip_mask

            beta_col = None
            for col in ct_bulk_decon_df.columns:
                if col.startswith("Beta") and col.endswith("_{}:GT".format(b_ct)):
                    beta_col = col
                    break

            print("\tPrepare plot df.")
            plot_df = pd.DataFrame({"sn_zscore": sn_df["OverallZScore"] * flip_mask,
                                    "b_zscore": ct_bulk_decon_df["OverallZScore"],
                                    "decon_beta": ct_bulk_decon_df[beta_col],
                                    "FDR": ct_bulk_decon_df["{}_FDR".format(b_ct)]})

            print("\tPlotting row 1.")
            self.plot(df=plot_df,
                      fig=fig,
                      ax=axes[0, col_index],
                      x="sn_zscore",
                      y="b_zscore",
                      xlabel="",
                      ylabel="bulk decon [z-score]",
                      name="{} VS {}".format(sn_ct, b_ct),
                      color=color)

            print("\tPlotting row 2.")
            self.plot(df=plot_df.loc[plot_df["FDR"] < 0.05, :],
                      fig=fig,
                      ax=axes[1, col_index],
                      x="sn_zscore",
                      y="b_zscore",
                      xlabel="",
                      ylabel="bulk decon [z-score]",
                      name="",
                      color=color)

            print("\tPlotting row 3.")
            self.plot(df=plot_df.loc[plot_df["FDR"] < 0.05, :],
                      fig=fig,
                      ax=axes[2, col_index],
                      x="sn_zscore",
                      y="decon_beta",
                      xlabel="single-nucleus [z-score]",
                      ylabel="bulk decon [log(abs(beta+1))]",
                      name="",
                      color=color)

        fig.savefig(os.path.join(self.outdir, "zscore_comparison.{}".format(self.extension)))
        plt.close()

    @staticmethod
    def log_modulus_beta(beta_df):
        df = beta_df.copy()
        data = []
        for _, row in beta_df.T.iterrows():
            values = []
            for beta in row:
                values.append(np.log(abs(beta)+1) * np.sign(beta))
            data.append(values)
        new_df = pd.DataFrame(data, index=df.columns, columns=df.index)

        return new_df.T

    @staticmethod
    def pvalue_to_zscore(pvalue_df):
        df = pvalue_df.copy()
        data = []
        for _, row in pvalue_df.T.iterrows():
            zscores = []
            for pvalue in row:
                if pvalue > (1.0 - 1e-16):
                    zscores.append(-8.209536151601387)
                elif pvalue < 1e-323:
                    zscores.append(-8.209536151601387)
                else:
                    zscores.append(stats.norm.isf(pvalue))
            data.append(zscores)
        zscore_df = pd.DataFrame(data, index=df.columns, columns=df.index)

        return zscore_df.T

    def create_all_eQTL_df(self, sn_df, bulk_decon_df, b_ct):
        overlap = set(sn_df.index).intersection(set(bulk_decon_df.index))

        df1 = sn_df.loc[overlap, :].copy()
        df2 = bulk_decon_df.loc[overlap, :].copy()

        flip_mask = self.create_flip_mask(df1, df2)

        df = pd.DataFrame({"x": df1["OverallZScore"] * flip_mask,
                           "y": df2[b_ct]})

        return df

    @staticmethod
    def create_flip_mask(df1, df2):
        if df1.shape[0] != df2.shape[0]:
            print("Unequal rows.")
            exit()

        flip_mask = []
        for i in range(len(df1.shape[0])):
            if df1.iloc[i, :]["SNPType"].split("/")[1] != df2.iloc[i, :]["AlleleAssessed"]:
                flip_mask.append(-1)
            else:
                flip_mask.append(1)

        return flip_mask

    def plot(self, df, fig, ax, x="x", y="y", xlabel="", ylabel="", name="",
             color="#000000"):
        sns.despine(fig=fig, ax=ax)

        #coef, p = stats.spearmanr(subset[x], subset[y])
        coef, p = stats.pearsonr(df[x], df[y])

        sns.regplot(x=x, y=y, data=df,
                    scatter_kws={'facecolors': "#808080",
                                 'edgecolors': "#808080"},
                    line_kws={"color": color},
                    ax=ax
                    )

        ax.axhline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)
        ax.axvline(0, ls='--', color="#D7191C", alpha=0.3, zorder=-1)

        ax.text(0.5, 1.1, name.replace("_", " "),
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02, "N = {}".format(df.shape[0]),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        handles = []
        handles.append(mpatches.Patch(color=color, label="{} [{:.2f}]".format(name.replace("_VS_", " / "), coef)))
        ax.legend(handles=handles)


if __name__ == '__main__':
    m = main()
    m.start()
