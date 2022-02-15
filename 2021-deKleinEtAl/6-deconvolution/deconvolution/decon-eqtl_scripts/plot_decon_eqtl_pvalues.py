#!/usr/bin/env python3

"""
File:         plot_decon_eqtl_pvalues.py
Created:      2022/02/15
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
import argparse
import math
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

"""
Syntax:
./plot_decon_eqtl_pvalues.py \
    -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron/deconvolutionResults.txt.gz \
    -o 2022-01-26-CortexEUR-cis-ForceNormalised-MAF5-4SD-CompleteConfigs-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron \
    -e png pdf
"""

# Metadata
__program__ = "Plot Decon-eQTL Nominal P-values"
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
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon_path = getattr(arguments, 'decon_path')
        self.outfile = getattr(arguments, 'outfile')
        self.extensions = getattr(arguments, 'extensions')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot_decon_eqtl_pvalues')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.palette = {
            "Astrocyte": "#D55E00",
            "EndothelialCell": "#CC79A7",
            "Excitatory": "#56B4E9",
            "Microglia": "#E69F00",
            "Oligodendrocyte": "#009E73",
            "OtherNeuron": "#0072B2"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-d",
                            "--decon_path",
                            type=str,
                            required=True,
                            help="The path to the deconvolution results matrix")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the output file")
        parser.add_argument("-e",
                            "--extensions",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        decon_df = self.load_file(self.decon_path, header=0, index_col=0)
        nominal_pvalues_df = decon_df.loc[:, [x for x in decon_df.columns if x.endswith("_pvalue")]]
        print(nominal_pvalues_df)
        nominal_pvalues_df.columns = [x.replace("_pvalue", "") for x in nominal_pvalues_df.columns]
        print(nominal_pvalues_df)
        cell_types = list(nominal_pvalues_df.columns)
        print("\tcell types: {}".format(", ".join(cell_types)))

        print("Plotting data")
        self.plot(df=nominal_pvalues_df,
                  cell_types=cell_types)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def plot(self, df, cell_types):
        nplots = len(cell_types)
        ncols = math.ceil(np.sqrt(nplots))
        nrows = math.ceil(nplots / ncols)

        sns.set(rc={'figure.figsize': (ncols * 8, nrows * 6)})
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharex='all',
                                 sharey='all')

        row_index = 0
        col_index = 0
        for i in range(ncols * nrows):
            if nrows == 1 and ncols == 1:
                ax = axes
            elif nrows == 1 and ncols > 1:
                ax = axes[col_index]
            elif nrows > 1 and ncols == 1:
                ax = axes[row_index]
            else:
                ax = axes[row_index, col_index]

            if i < len(cell_types):
                ct = cell_types[i]

                sns.despine(fig=fig, ax=ax)

                sns.histplot(data=df[ct],
                             ax=ax,
                             kde=True,
                             kde_kws={"cut": 0},
                             color=self.palette[ct])

                ax.set_title(ct,
                             fontsize=18,
                             fontweight='bold')
                ax.set_ylabel("count",
                              fontsize=14,
                              fontweight='bold')
                ax.set_xlabel("nominal p-value",
                              fontsize=14,
                              fontweight='bold')

                ax.annotate(
                    'N = {:,}'.format(df.shape[0]),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color=self.palette[ct],
                    fontsize=14,
                    fontweight='bold'
                )

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        # Add the main title.
        fig.suptitle("Decon-eQTL nominal p-values",
                     fontsize=25,
                     color="#000000",
                     weight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}.{}".format(self.outfile, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Decon path: {}".format(self.decon_path))
        print("  > Output filename: {}".format(self.outfile))
        print("  > Extensions: {}".format(", ".join(self.extensions)))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
