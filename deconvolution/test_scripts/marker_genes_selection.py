#!/usr/bin/env python3

"""
File:         marker_genes_selection.py
Created:      2020/04/08
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
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Marker Gene Selection"
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
        self.excel_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/41598_2018_27293_MOESM3_ESM.xlsx"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the excel file")
        df = pd.read_excel(self.excel_path, header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.excel_path),
                                      df.shape))

        # Subset.
        df = df.loc[df["Cell_type_mentions"] > 0, :]

        # Remap.
        df["Celltype"] = df["Celltype"].map({"opc": "oligodendrocyte precursor cell",
                                             "oli": "oligodendrocyte",
                                             "neu": "neuron",
                                             "mic": "microglia",
                                             "end": "endothelial cell",
                                             "ast": "astrocyte"})

        # Order.
        df.sort_values(["Celltype", "Cell_type_mentions"], inplace=True,
                       ascending=False)
        df['x'] = df.groupby('Celltype').cumcount()

        # Plot.
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")

        g = sns.lineplot(x="x", y="Cell_type_mentions", hue="Celltype", data=df)
        ax.axvline(5, ls='--', color="#000000", zorder=-1)

        g.set_title("McKenzie et al. 2018 - Marker Gene Mentions",
                    fontsize=12,
                    fontweight='bold')
        g.set_ylabel("cell type mentions",
                     fontsize=10,
                     fontweight='bold')
        g.set_xlabel("")
        plt.legend(loc=7, prop={'size': 5})
        plt.setp(ax.get_legend().get_texts(), fontsize='10')
        plt.setp(ax.get_legend().get_title(), fontsize='12')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "mckenzie_cell_type_mentions.png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
