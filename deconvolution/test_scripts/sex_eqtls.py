#!/usr/bin/env python3

"""
File:         sex_eqtls.py
Created:      2020/05/01
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
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Sex eQTLs"
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
        # Interaction matrix files.
        self.inter_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/custom_interaction_analyser/output_backup_single_core_jobs/sex_interaction_table_with_indices.txt.gz"
        self.eqtl_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/matrix_preparation/output/combine_eqtlprobes/eQTLprobes_combined.txt.gz"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading dataframes.")
        inter_df = pd.read_csv(self.inter_path, sep="\t", header=0, index_col=0)
        print("\t{}".format(inter_df.shape))
        eqtl_df = pd.read_csv(self.eqtl_path, sep="\t", header=0)
        print("\t{}".format(eqtl_df.shape))

        inter_df = inter_df.T

        inter_df.reset_index(inplace=True)
        inter_df["index"] = inter_df["index"].astype(str)
        eqtl_df.reset_index(inplace=True)
        eqtl_df["index"] = eqtl_df["index"].astype(str)

        df = inter_df.merge(eqtl_df, how='inner', left_on=['index', '-'], right_on=['index', 'SNPName'])
        df["SEX"] = df["SEX"].astype(float)
        sex_df = df.loc[df.loc[:, "SEX"] > 2.89, :].copy()
        sex_df.sort_values(by="SEX", inplace=True)
        print(sex_df)
        genes = list(sex_df["HGNCName"].values)
        print(genes)

if __name__ == '__main__':
    m = main()
    m.start()
