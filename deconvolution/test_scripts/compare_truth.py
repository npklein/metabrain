#!/usr/bin/env python3

"""
File:         gradient_plot.py
Created:      2020/04/15
Last Changed: 2020/06/02
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
__program__ = "Gradient Plot"
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
       self.sc_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/single_cell_counts.txt.gz"
       self.ihc_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/test_scripts/AMP-AD/IHC_counts.txt.gz"
       self.outdir = str(Path(__file__).parent.parent)
       self.palette = {
           "Neuron": "#0072B2",
           "Oligodendrocyte": "#009E73",
           "EndothelialCell": "#CC79A7",
           "Microglia": "#E69F00",
           "Macrophage": "#E69F00",
           "Astrocyte": "#D55E00",
           "Pericytes": "#808080"
       }

    def start(self):
        print("Load the single-cell file.")
        sc_df = pd.read_csv(self.sc_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.sc_path),
                                      sc_df.shape))

        print("Load the IHC file.")
        ihc_df = pd.read_csv(self.ihc_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.ihc_path),
                                      ihc_df.shape))

        overlap = np.intersect1d(sc_df.index, ihc_df.index)
        sc_df = sc_df.loc[overlap, :]
        ihc_df = ihc_df.loc[overlap, :]

        sc_df.index.name = "-"
        ihc_df.index.name = "-"

        if not sc_df.index.equals(ihc_df.index):
            print("Invalid order")
            exit()

        # Create the comparison dataframe.
        cell_types = np.intersect1d(sc_df.columns, ihc_df.columns)
        data = []
        for i, sample in enumerate(overlap):
            sc = sc_df.iloc[i, :]
            ihc = ihc_df.iloc[i, :]
            for cell in cell_types:
                data.append([sample, sc[cell], ihc[cell], cell])
        df = pd.DataFrame(data, columns=['-', 'x', 'y', 'hue'])
        df.set_index('-', inplace=True)

        # Plot.
        df['color'] = df['hue'].map(self.palette)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        groups = df['hue'].unique()

        coefs = {}
        for group in groups:
            subset = df.loc[df['hue'] == group, :].copy()
            color = self.palette[group]

            coef, p = stats.spearmanr(subset["x"], subset["y"])
            coefs[group] = "r = {:.2f}".format(coef)

            sns.regplot(x="x", y="y", data=subset,
                        scatter_kws={'facecolors': subset['color'],
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.plot([-0.1, 1.1], [-0.1, 1.1], ls="--", c=".3")

        ax.text(0.5, 1.1, "single-cell versus IHC",
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02, "N = {}".format(len(overlap)),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel("IHC counts",
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel("single-cell counts",
                      fontsize=14,
                      fontweight='bold')

        handles = []
        for key, color in self.palette.items():
            if key in groups and key in coefs.keys():
                handles.append(mpatches.Patch(color=color, label="{} [{}]".format(key, coefs[key])))
        ax.legend(handles=handles)

        fig.savefig(os.path.join(self.outdir, "ground_truth_comparison.png"))
        plt.close()

if __name__ == '__main__':
    m = main()
    m.start()
