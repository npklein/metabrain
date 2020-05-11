"""
File:         deconvolution_covariate_comparison.py
Created:      2020/04/15
Last Changed: 2020/04/28
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
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir, p_value_to_symbol


class DeconvolutionCovariateComparison:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'deconvolution_covariate_comparison')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.cov_df = dataset.get_cov_df()

    def start(self):
        print("Plotting deconvolution convariate comparison.")
        self.print_arguments()
        corr_df, pval_df = self.correlate(self.cov_df)
        self.plot(corr_df, pval_df, self.outdir)

    @staticmethod
    def correlate(df):
        print("Creating the correlation matrix.")
        corr_df = pd.DataFrame(np.nan, index=df.index, columns=df.index)
        pval_df = pd.DataFrame("", index=df.index, columns=df.index)
        for i, row1 in enumerate(df.index):
            for j, row2 in enumerate(df.index):
                if i >= j:
                    coef, p = stats.spearmanr(df.loc[row1, :], df.loc[row2, :])
                    corr_df.loc[row1, row2] = coef
                    pval_df.loc[row1, row2] = p_value_to_symbol(p)

        return corr_df, pval_df

    @staticmethod
    def plot(corr_df, pval_df, outdir):
        print("Plotting")

        indices = [(93, 98, "McKenzie\nMG", "McKenzie_"),
                   (98, 103, "McKenzie\nMG", "McKenzie_"),
                   (103, 108, "McKenzie\nMG", "McKenzie_"),
                   (108, 113, "McKenzie\nMG", "McKenzie_"),
                   (113, 118, "McKenzie\nMG", "McKenzie_"),
                   (118, 123, "CellMap\nPCA", "CellMapPCA_"),
                   (123, 128, "CellMap\nNMF", "CellMapNMF_"),
                   (128, 133, "CellMap\nNNLS", "CellMapNNLS_")]

        cmap = plt.cm.RdBu_r
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        heatmapkws = dict(square=False, cbar=False, cmap=cmap, fmt='',
                          linewidths=1.0, center=0, vmin=-1, vmax=1,
                          annot_kws={"size": 12, "color": "#808080"})

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=len(indices), nrows=len(indices),
                                 figsize=(29, 21))
        plt.subplots_adjust(left=0.2, right=0.87, bottom=0.2, top=0.95,
                            wspace=0.1, hspace=0.1)

        for i in range(len(indices)):
            (ra, rb, ylabel, yremove) = indices[i]
            for j in range(len(indices)):
                (ca, cb, xlabel, xremove) = indices[j]
                ax = axes[i, j]

                if i >= j:
                    print("\tPlotting axes[{}, {}]".format(i, j))
                    xticklabels = False
                    if i == len(indices) - 1:
                        xticklabels = True
                    yticklabels = False
                    if j == 0:
                        yticklabels = True

                    data_subset = corr_df.iloc[ra:rb, ca:cb]
                    anmnot_subset = pval_df.iloc[ra:rb, ca:cb]

                    sns.heatmap(data_subset,
                                annot=anmnot_subset,
                                xticklabels=xticklabels,
                                yticklabels=yticklabels,
                                ax=ax,
                                **heatmapkws)

                    if xticklabels:
                        new_xticks = [x.replace(xremove, '').replace("_", " ") for x in data_subset.columns]
                        ax.set_xticklabels(new_xticks, fontsize=14, rotation=90)
                        ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')

                    if yticklabels:
                        new_yticks = [y.replace(yremove, '').replace("_", " ") for y in data_subset.index]
                        ax.set_yticklabels(new_yticks, fontsize=14, rotation=0)
                        ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
                else:
                    ax.set_axis_off()

        cax = fig.add_axes([0.9, 0.2, 0.01, 0.7])
        sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.align_ylabels(axes[:, 0])
        fig.align_xlabels(axes[len(indices) - 1, :])
        fig.colorbar(sm, cax=cax)
        fig.suptitle('Deconvolution Covariate Correlations', fontsize=40, fontweight='bold')
        fig.savefig(os.path.join(outdir, "deconvolution_covariate_comparison.pdf"), format='pdf', dpi=600)
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
