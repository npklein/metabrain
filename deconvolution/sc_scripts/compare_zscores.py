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
        self.sn_infolder = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-03-ROSMAP-scRNAseq/"
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

        print("Preprocessing deconvolution data frame")
        decon_work_df = self.filter_df(decon_df, prefix="Beta", suffix=":GT")
        decon_work_df.columns = [x.split("_")[1] for x in decon_work_df.columns]
        decon_work_df = self.log_modulus_beta(decon_work_df)
        ylabel_suffix = "log modulus beta"

        # decon_work_df = self.filter_df(decon_df, suffix="_FDR")
        # ylabel_suffix = "OverallZScore"

        # decon_work_df = self.filter_df(decon_df, suffix="_pvalue")
        # decon_work_df = self.pvalue_to_zscore(decon_work_df)
        # ylabel_suffix = "z_score"

        probe_names = []
        snp_names = []
        for index in decon_work_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        decon_work_df["ProbeName"] = probe_names
        decon_work_df["SNPName"] = snp_names
        decon_work_df.reset_index(drop=True, inplace=True)

        print("Merging results")
        bulk_decon_df = bulk_df.merge(decon_work_df,
                                      left_on=["SNPName", "ProbeName"],
                                      right_on=["SNPName", "ProbeName"])
        del bulk_df

        bulk_decon_df.index = bulk_decon_df["SNPName"] + "_" + bulk_decon_df["ProbeName"]
        print(bulk_decon_df)

        for sn_ct, b_ct, color in self.cell_types:
            print("Analyzing {} / {}".format(sn_ct, b_ct))

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

            print("\tPrepare plot df.")
            plot_df = pd.DataFrame({"x": sn_df["OverallZScore"] * flip_mask,
                                    "y": ct_bulk_decon_df[b_ct]})

            print("\tPlotting.")
            self.plot(plot_df,
                      xlabel="single-nucleus [OverallZScore]",
                      ylabel="bulk deconvolution [{}]".format(ylabel_suffix),
                      name="{}_VS_{}".format(sn_ct, b_ct),
                      color=color)

    @staticmethod
    def filter_df(df, prefix=None, suffix=None):
        cols = []
        for col in df.columns:
            if (prefix is not None and suffix is None and col.startswith(prefix)) or \
                    (prefix is None and suffix is not None and col.endswith(suffix)) or \
                    (prefix is not None and suffix is not None and col.startswith(prefix) and col.endswith(suffix)):
                cols.append(col)
        if len(cols) <= 0:
            print("Error, no columns selected.")
            exit()
        df = df[cols]
        for string in [prefix, suffix]:
            if string is not None:
                df.columns = [x.replace(string, "") for x in cols]
        return df

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

    def plot(self, df, xlabel="", ylabel="", name="", color="#000000"):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        #coef, p = stats.spearmanr(subset["x"], subset["y"])
        coef, p = stats.pearsonr(df["x"], df["y"])

        sns.regplot(x="x", y="y", data=df,
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

        fig.savefig(os.path.join(self.outdir, "{}_regression.{}".format(name, self.extension)))
        plt.close()









if __name__ == '__main__':
    m = main()
    m.start()
