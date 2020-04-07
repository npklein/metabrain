#!/usr/bin/env python3

"""
File:         cellmap_pc_correlation.py
Created:      2020/04/07
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
        self.profile_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/perform_celltype_pca/celltype_pcs.txt.gz"
        self.marker_genes = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-deconvolution/matrix_preparation/output/create_matrices/marker_genes.txt.gz"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the CellMap profile PC1")
        df1 = pd.read_csv(self.profile_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      df1.shape))

        print("Load the marker genes")
        df2 = pd.read_csv(self.marker_genes, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.marker_genes),
                                      df2.shape))
        df2.drop_duplicates(inplace=True)
        df2.sort_index(inplace=True)

        print("Combining dataframes")
        df = pd.concat([df1, df2], axis=0)
        print("\tNew shape: {}".format(df.shape))

        corr_df = pd.DataFrame(np.nan, index=df.index, columns=df.index)
        pval_df = pd.DataFrame("", index=df.index, columns=df.index)
        for i, row1 in enumerate(df.index):
            for j, row2 in enumerate(df.index):
                if i >= j:
                    coef, p = stats.spearmanr(df.loc[row1, :], df.loc[row2, :])
                    corr_df.loc[row1, row2] = coef
                    pval_df.loc[row1, row2] = self.p_value_to_symbol(p)

        sns.set(color_codes=True)
        g = sns.clustermap(corr_df, center=0, cmap="RdBu_r",
                           row_cluster=False, col_cluster=False,
                           yticklabels=True, xticklabels=True,
                           vmin=-1, vmax=1, annot=pval_df, fmt='',
                           annot_kws={"size": 6, "color": "#808080"},
                           figsize=(12, 9))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10, rotation=0))
        plt.setp(
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(),
                                         fontsize=10, rotation=90))
        g.fig.suptitle('CellMap Profile PC1 Spearman Correlation', fontsize=22)
        g.savefig(os.path.join(self.outdir, "cellmap_profile_pc1_corr.png"))
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
