"""
File:         deconvolution_covariate_comparison.py
Created:      2020/04/15
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
        self.methods = dataset.get_methods()

    def start(self):
        print("Plotting deconvolution convariate comparison.")
        self.print_arguments()
        df = self.filter()
        corr_df, pval_df = self.correlate(df)
        self.plot(corr_df, pval_df, self.outdir)

    def filter(self):
        print("Filtering the covariate dataframe.")
        subset = []
        for index in self.cov_df.index:
            for method in self.methods:
                if method in index:
                    subset.append(index)
                    break
        df = self.cov_df.copy()
        df = df.loc[subset, :]
        print("\tShape after subsetting: {}".format(df.shape))
        return df

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

        indices = [(0, 5, "McKenzie", "McKenzie_"),
                   (5, 10, "McKenzie", "McKenzie_"),
                   (10, 15, "McKenzie", "McKenzie_"),
                   (15, 20, "McKenzie", "McKenzie_"),
                   (20, 25, "McKenzie", "McKenzie_"),
                   (25, 30, "CellMap", "CellMap_"),
                   (30, 35, "NNLS", "NNLS_")]

        cmap = plt.cm.RdBu_r
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        heatmapkws = dict(square=False, cbar=False, cmap=cmap, fmt='',
                          linewidths=1.0, center=0, vmin=-1, vmax=1,
                          annot_kws={"size": 12, "color": "#808080"})

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=7, nrows=7, figsize=(25, 18))
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
                    if i == 6:
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
                        new_xticks = [x.replace(xremove, '') for x in data_subset.columns]
                        ax.set_xticklabels(new_xticks, fontsize=14, rotation=90)
                        ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')

                    if yticklabels:
                        new_yticks = [y.replace(yremove, '') for y in data_subset.index]
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
        fig.savefig(os.path.join(outdir, "deconvolution_covariate_comparison.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Methods: {}".format(self.methods))
        print("  > Output directory: {}".format(self.outdir))
        print("")