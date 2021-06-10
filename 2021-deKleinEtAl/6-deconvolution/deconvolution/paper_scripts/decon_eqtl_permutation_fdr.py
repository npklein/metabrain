#!/usr/bin/env python3

"""
File:         decon_eqtl_permutation_fdr.py
Created:      2021/06/07
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
from bisect import bisect_left
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

"""
Syntax:
./decon-eqtl_permutation_fdr.py
"""

# Metadata
__program__ = "Decon-eQTL Permutation FDR"
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
        self.workdir = os.path.join(str(Path(__file__).parent.parent.absolute()), "decon_eqtl_permutation_per_cohort")

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self):
        self.print_arguments()

        # Load the real data.
        print("Loading real results.")
        real_path = os.path.join(self.workdir, "decon_output", "real", "deconvolutionResults.csv")
        real_df = self.load_file(real_path, header=0, index_col=0)
        index = real_df.index.tolist()
        columns = [x for x in real_df.columns if "pvalue" in x]
        real_df = real_df.loc[:, columns]
        print(real_df)
        real_m = real_df.to_numpy()
        del real_df

        print("Plotting real p-values.")
        self.histplot(data=pd.DataFrame(real_m.flatten(), columns=["p-values"]),
                      x="p-values", filename="real_pvalues")

        # Load the permutated data.
        print("Loading permuted genotype results.")
        perm_folders = [os.path.basename(x) for x in glob.glob(os.path.join(self.workdir, "decon_output", "*"))]
        perm_folders.remove('real')
        n_permutations = len(perm_folders)
        perm_m = np.empty((real_m.shape[0], real_m.shape[1], n_permutations), dtype=np.float64)
        for i, perm_folder in enumerate(perm_folders):
            perm_path = os.path.join(self.workdir, "decon_output", perm_folder, "deconvolutionResults.csv")
            perm_df = self.load_file(perm_path, header=0, index_col=0)
            # TODO match indices
            perm_m[:, :, i] = perm_df.loc[:, columns].to_numpy()
            del perm_df

        print("Plotting real p-values.")
        self.histplot(data=pd.DataFrame(perm_m.flatten(), columns=["perm. p-values"]),
                      x="perm. p-values", filename="perm_pvalues")

        # Calculate the fdr.
        print("Calculating permutation FDR per ieQTL.")
        ieqtl_fdr_m = np.empty_like(real_m, dtype=np.float64)
        for i in range(real_m.shape[0]):
            for j in range(real_m.shape[1]):
                ieqtl_fdr_m[i, j] = np.sum(perm_m[i, j, :] < real_m[i, j]) / n_permutations

        ieqtl_fdr_df = pd.DataFrame(ieqtl_fdr_m, index=index, columns=[x.replace("pvalue", "FDR") for x in columns])
        print(ieqtl_fdr_df)
        ieqtl_fdr_path = os.path.join(self.workdir, "deconvolutionResults_{}perm_perIeQTLFDR.txt.gz".format(n_permutations))
        self.save_file(ieqtl_fdr_df, outpath=ieqtl_fdr_path)

        print("Calculating permutation FDR over all permutations.")
        comb_fdr_m = np.empty_like(real_m, dtype=np.float64)
        for i in range(real_m.shape[0]):
            for j in range(real_m.shape[1]):
                comb_fdr_m[i, j] = np.sum(perm_m < real_m[i, j]) / np.size(perm_m)
                # sorteer echte pwaardes en gepermuteerde p waardes
                # tel de fracties van p waardes < echte p waarde in de null disitrubtie (count / n permutations)
                # ene proportie en de andere

        comb_fdr_df = pd.DataFrame(comb_fdr_m, index=index, columns=[x.replace("pvalue", "FDR") for x in columns])
        print(comb_fdr_df)
        comb_fdr_path = os.path.join(self.workdir, "deconvolutionResults_{}perm_CombinedPermFDR.txt.gz".format(n_permutations))
        self.save_file(comb_fdr_df, outpath=comb_fdr_path)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def histplot(self, data, x="x", xlabel=None, title="", filename="plot"):
        if xlabel is None:
            xlabel = x

        sns.set(rc={'figure.figsize': (10, 7.5)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        # plot in sections.
        sns.histplot(data=data, x=x, ax=ax)

        ax.set_title(title,
                     fontsize=18,
                     fontweight='bold')
        ax.set_ylabel("count",
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        fig.savefig(os.path.join(self.outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Work directory: {}".format(self.workdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
