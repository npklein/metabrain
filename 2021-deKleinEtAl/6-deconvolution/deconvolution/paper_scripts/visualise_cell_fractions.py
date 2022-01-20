#!/usr/bin/env python3

"""
File:         visualise_cell_fractions.py
Created:      2020/11/24
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

# Metadata
__program__ = "Visualise Cell Fractions"
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

"""
Syntax:
./visualise_cell_fractions.py -cf ../matrix_preparation/CortexEUR-cis-Uncorrected/perform_deconvolution/deconvolution_table.txt.gz -e png pdf -o CortexEUR-cis-Uncorrected-PsychENCODEProfile

./visualise_cell_fractions.py -cf ../matrix_preparation/CortexEUR-cis-Uncorrected-NoENA-NoMDSOutlier/perform_deconvolution/deconvolution_table.txt.gz -e png pdf -o CortexEUR-cis-Uncorrected-NoENA-NoMDSOutlier-PsychENCODEProfile

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz -e png pdf -o 2021-12-07-CortexEUR-PsychENCODEProfileSummedNoDev

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_NoInhibitory.txt.gz -e png pdf -o 2021-12-07-CortexEUR-PsychENCODEProfileSummedNoDevNoInhibitory

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron.txt.gz -e png pdf -o 2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved-deconvolutionTable-InhibitorySummedWithOtherNeuron

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table.txt.gz -e png pdf -o 2021-12-22-CortexAFR-replicationOfCortexEUR20211207-cis-ProbesWithZeroVarianceRemove-InhibitorySummedWithOtherNeuron

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/perform_deconvolution/deconvolution_table_InhibitorySummedWithOtherNeuron_DatasetNormalised.txt.gz -e png pdf -o 2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved-deconvolutionTable-InhibitorySummedWithOtherNeuron-NNLSBetasDatasetNormnalised

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-19-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -e png pdf -o 2022-01-19-CortexEUR-cis-deconvolutionTable

./visualise_cell_fractions.py -cf /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2022-01-20-CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -e png pdf -o 2022-01-20-CortexEUR-cis-deconvolutionTable
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.cf_path = getattr(arguments, 'cell_fractions')
        self.extensions = getattr(arguments, 'extension')
        self.output_filename = getattr(arguments, 'output')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'visualise_cell_fractions')
        self.palette = {
            'Adult-Ex1': '#56B4E9',
            'Adult-Ex2': '#56B4E9',
            'Adult-Ex3': '#56B4E9',
            'Adult-Ex4': '#56B4E9',
            'Adult-Ex5': '#56B4E9',
            'Adult-Ex6': '#56B4E9',
            'Adult-Ex7': '#56B4E9',
            'Adult-Ex8': '#56B4E9',
            'Adult-In1': '#0072B2',
            'Adult-In2': '#0072B2',
            'Adult-In3': '#0072B2',
            'Adult-In4': '#0072B2',
            'Adult-In5': '#0072B2',
            'Adult-In6': '#0072B2',
            'Adult-In7': '#0072B2',
            'Adult-In8': '#0072B2',
            'Adult-Microglia': '#E69F00',
            'Adult-OPC': '#1b8569',
            'Adult-Endothelial': '#CC79A7',
            'Adult-Astrocytes': '#D55E00',
            'Adult-Oligo': '#009E73',
            'Adult-OtherNeuron': '#2690ce',
            'Dev-Replicating': '#000000',
            'Dev-Quiescent': '#808080',
            "Excitatory": "#56B4E9",
            "Inhibitory": "#2690ce",
            'OtherNeuron': '#0072B2',
            "Oligodendrocyte": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Astrocyte": "#D55E00"
        }

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # Set the right pdf font for exporting.
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-cf",
                            "--cell_fractions",
                            type=str,
                            required=True,
                            help="The path to the cell fractions matrix.")
        parser.add_argument("-e",
                            "--extension",
                            nargs="+",
                            type=str,
                            choices=["png", "pdf", "eps"],
                            default=["png"],
                            help="The figure file extension. "
                                 "Default: 'png'.")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            required=True,
                            help="The name of the output image.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data.")
        cf_df = self.load_file(self.cf_path)
        order = list(cf_df.columns)

        print("Preprocessing data.")
        cc_dfm = cf_df.melt()
        print(cc_dfm)

        print("Plotting.")
        self.create_boxplot(df=cc_dfm,
                            hue="variable",
                            order=order,
                            palette=self.palette,
                            ylabel="cell fraction %",
                            name=self.output_filename)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  low_memory=True):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows, low_memory=low_memory)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    def create_boxplot(self, df, x="variable", y="value", hue=None, xlabel="",
                       ylabel="", name="", order=None, palette=None):
        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        sns.violinplot(x=x,
                       y=y,
                       hue=hue,
                       data=df,
                       order=order,
                       palette=palette,
                       cut=0,
                       dodge=False,
                       ax=ax)

        plt.setp(ax.collections, alpha=.75)

        sns.boxplot(x=x,
                    y=y,
                    hue=hue,
                    data=df,
                    order=order,
                    whis=np.inf,
                    color="white",
                    dodge=False,
                    ax=ax)

        if ax.get_legend() is not None:
            ax.get_legend().remove()

        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')

        for extension in self.extensions:
            fig.savefig(os.path.join(self.outdir, "{}_boxplot.{}".format(name, extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Cell fractions path: {}".format(self.cf_path))
        print("  > Output filename: {}".format(self.output_filename))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
