#!/usr/bin/env python3

"""
File:         deconvolution_covariate_comparison.py
Created:      2020/04/07
Last Changed: 2020/04/14
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
from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "CellMap Profiles"
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


class main():
    def __init__(self):
        self.covariates_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/output/create_cov_matrix/covariates_table.txt.gz"
        self.methods = ["McKenzie", "CellMap", "NNLS"]
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the covariate matrix.")
        df = pd.read_csv(self.covariates_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.covariates_path),
                                      df.shape))

        subset = []
        for index in df.index:
            for method in self.methods:
                if method in index:
                    subset.append(index)
                    break
        df = df.loc[subset, :]
        print("Shape after subsetting: {}".format(df.shape))

        corr_df = pd.DataFrame(np.nan, index=df.index, columns=df.index)
        pval_df = pd.DataFrame("", index=df.index, columns=df.index)
        for i, row1 in enumerate(df.index):
            for j, row2 in enumerate(df.index):
                if i >= j:
                    coef, p = stats.spearmanr(df.loc[row1, :], df.loc[row2, :])
                    corr_df.loc[row1, row2] = coef
                    pval_df.loc[row1, row2] = self.p_value_to_symbol(p)

        # Plot.
        left = 0.15
        right = 0.87
        bottom = 0.15
        top = 0.9
        cmap = plt.cm.RdBu_r
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        heatmapkws = dict(square=False, cbar=False, cmap=cmap, fmt='',
                          linewidths=1.0, center=0, vmin=-1, vmax=1,
                          annot_kws={"size": 12, "color": "#808080"})

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=7, nrows=7, figsize=(25, 18))
        plt.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                            wspace=0.1, hspace=0.1)

        indices = list(zip([x for x in range(0, 35, 5)], [x for x in range(5, 40, 5)]))
        for i in range(7):
            ri = indices[i]
            for j in range(7):
                ci = indices[j]
                ax = axes[i, j]

                if i >= j:
                    xticklabels = False
                    if i == 6:
                        xticklabels = True
                    yticklabels = False
                    if j == 0:
                        yticklabels = True

                    data_subset = corr_df.iloc[ri[0]:ri[1], ci[0]:ci[1]]
                    anmnot_subset = pval_df.iloc[ri[0]:ri[1], ci[0]:ci[1]]

                    sns.heatmap(data_subset,
                                annot=anmnot_subset,
                                xticklabels=xticklabels,
                                yticklabels=yticklabels,
                                ax=ax,
                                **heatmapkws)

                    if xticklabels:
                        new_labels = [' '.join(x.split("_")[1:]) for x in data_subset.columns]
                        ax.set_xticklabels(new_labels, fontsize=14, rotation=90)
                        ax.set(xlabel=corr_df.index[ci[0]].split("_")[0])

                    if yticklabels:
                        new_labels = [' '.join(str(x).split("_")[1:]) for x in data_subset.index]
                        ax.set_yticklabels(new_labels, fontsize=14, rotation=0)
                        ax.set(ylabel=corr_df.index[ri[0]].split("_")[0])
                else:
                    ax.set_axis_off()

        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, cax=cax)

        # sns.set(color_codes=True)
        # g = sns.clustermap(corr_df, center=0, cmap="RdBu_r",
        #                    row_cluster=False, col_cluster=False,
        #                    row_colors=colors, col_colors=colors,
        #                    yticklabels=True, xticklabels=True,
        #                    vmin=-1, vmax=1, annot=pval_df, fmt='',
        #                    annot_kws={"size": 6, "color": "#808080"},
        #                    figsize=(12, 9))
        # plt.setp(
        #     g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
        #                                  fontsize=10, rotation=0))
        # plt.setp(
        #     g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),
        #                                  fontsize=10, rotation=90))
        fig.suptitle('Deconvolution Covariate Correlations', fontsize=40)
        fig.savefig(os.path.join(self.outdir, "deconvolution_covariate_comparison.png"))
        plt.close()

    @staticmethod
    def p_value_to_symbol(p_value):
        output = ""
        try:
            if p_value > 0.05:
                output = "NS"
            elif 0.05 >= p_value > 0.01:
                output = "*"
            elif 0.01 >= p_value > 0.001:
                output = "**"
            elif 0.001 >= p_value > 0.0001:
                output = "***"
            elif p_value <= 0.0001:
                output = "****"
        except TypeError:
            pass

        return output


if __name__ == '__main__':
    m = main()
    m.start()
