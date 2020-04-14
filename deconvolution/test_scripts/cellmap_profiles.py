#!/usr/bin/env python3

"""
File:         cellmap_profiles.py
Created:      2020/04/06
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
import scipy.stats as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

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
        self.profile_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/CellMap_brain_celltype_avgCPM.txt"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Load the profile")
        df = pd.read_csv(self.profile_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.profile_path),
                                      df.shape))

        print("Filtering genes with 0 std")
        df = df.loc[df.std(axis=1) > 0, :]
        df = df.T
        normalized_df = (df-df.mean()) / df.std()
        print("\tNew shape: {}".format(normalized_df.shape))


        sns.set(color_codes=True)
        g = sns.clustermap(normalized_df, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=False,
                           figsize=(12, 9))
        plt.setp(g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.fig.suptitle('CellMap Profile')
        plt.tight_layout()
        g.savefig(os.path.join(self.outdir, "cellmap_profile.png"))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
