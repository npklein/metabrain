#!/usr/bin/env python3

"""
File:         julienbryois2021_replication_plot3.py
Created:      2022/01/24
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
from statsmodels.stats import multitest
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.

"""
Syntax:
./julienbryois2021_replication_plot3.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-19-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-ExprMatrixCovariatesRemovedOLS-NegativeToZero/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL -o 2022-01-19-CortexEUR-cis-NormalisedMAF5-LimitedConfigs-PsychENCODEProfile-NoDev-InhibitorySummedWithOtherNeuron-ExprMatrixCovariatesRemovedOLS-NegativeToZero

./julienbryois2021_replication_plot3.py -dd /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron-GT5SD/deconvolutionResults.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/2021-12-07-CortexEUR-cis-ProbesWithZeroVarianceRemoved/create_matrices/genotype_alleles.txt.gz -dn MetaBrain_Decon-eQTL -o 2022-01-21-CortexEUR-cis-NegativeToZero-DatasetAndRAMCorrected-InhibitorySummedWithOtherNeuron-GT5SD

"""

# Metadata
__program__ = "Julien Bryois et al. 2021 Replication Plot 3"
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
        self.discovery_data_path = getattr(arguments, 'discovery_data_path')
        self.discovery_alleles_path = getattr(arguments, 'discovery_alleles_path')
        self.discovery_name = getattr(arguments, 'discovery_name')
        self.out_filename = getattr(arguments, 'outfile')
        self.replication_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/julienbryois2021"
        self.replication_name = "JulienBryois2021"

        # Set variables.
        outdir = os.path.join(str(Path(__file__).parent.parent), 'julienbryois2021_replication_plot')
        self.outdir = os.path.join(outdir, 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.colormap = {
            "Excitatory": "#56B4E9",
            "ExcitatoryNeurons": "#56B4E9",
            "Inhibitory": "#0072B2",
            "InhibitoryNeurons": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "Oligodendrocytes": "#009E73",
            "OPC": "#009E73",
            "OPCs...COPs": "#009E73",
            "OPCsCOPs": "#009E73",
            "EndothelialCell": "#CC79A7",
            "EndothelialCells": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Astrocytes": "#D55E00",
            "Pericytes": "#808080",
            "OtherNeuron": "#0072B2"
        }

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
                            help="show program's version number and exit")
        parser.add_argument("-dd",
                            "--discovery_data_path",
                            type=str,
                            required=True,
                            help="The path to the discovery deconvolution "
                                 "results matrix")
        parser.add_argument("-da",
                            "--discovery_alleles_path",
                            type=str,
                            required=True,
                            help="The path to the discovery alleles matrix")
        parser.add_argument("-dn",
                            "--discovery_name",
                            type=str,
                            default="one",
                            help="The name for the discovery data")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=True,
                            help="The name of the outfile.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading data")
        discovery_df = self.load_file(self.discovery_data_path, header=0, index_col=None)
        discovery_df["SNPName"] = ["_".join(x.split("_")[1:]) for x in discovery_df.iloc[:, 0]]
        discovery_cell_types = [x.replace("_pvalue", "") for x in discovery_df if x.endswith("_pvalue")]

        discovery_alleles_df = self.load_file(self.discovery_alleles_path, header=0, index_col=None)
        discovery_alleles_df.columns = ["SNPName", "Alleles", "MinorAllele"]
        discovery_alleles_df["DeconAllele"] = [x.split("/")[1] for x in discovery_alleles_df["Alleles"]]
        snp_to_da_dict = dict(zip(discovery_alleles_df["SNPName"], discovery_alleles_df["DeconAllele"]))

        discovery_df["DeconAllele"] = discovery_df["SNPName"].map(snp_to_da_dict)
        discovery_df.index = ["{}_{}".format(x.split(".")[0], x.split(":")[2]) for x in discovery_df.iloc[:, 0]]
        discovery_df = discovery_df.iloc[:, 1:]
        discovery_df.reset_index(drop=False, inplace=True)

        discovery_fdr_df = discovery_df.loc[:, [x for x in discovery_df.columns if x.endswith("_FDR")]]
        discovery_fdr_df = (discovery_fdr_df < 0.05).astype(int)
        discovery_unique_ieqtls_df = discovery_fdr_df.loc[discovery_fdr_df.sum(axis=1) == 1, :]

        replication_df = self.load_file(os.path.join(self.replication_path, "JulienBryois2021SummaryStats.txt.gz"), header=0, index_col=0, nrows=None)
        replication_cell_types = [x.replace(" p-value", "") for x in replication_df if x.endswith(" p-value")]
        replication_df.reset_index(drop=False, inplace=True)

        indices = []
        replication_data = []
        annot_data = []
        for discovery_ct in discovery_cell_types:
            discovery_fdr_column = "{}_FDR".format(discovery_ct)
            discovery_beta_column = None
            for col in discovery_df:
                if discovery_ct in col and col.endswith(":GT"):
                    discovery_beta_column = col
                    break

            # Select the discovery ieQTLs.
            ieqtls = list(discovery_unique_ieqtls_df.loc[discovery_unique_ieqtls_df[discovery_fdr_column] == 1, :].index)
            # ieqtls = list(discovery_df.loc[discovery_df[discovery_fdr_column] < 0.05, :].index)
            discovery_subset_df = discovery_df.loc[ieqtls, ["index", discovery_beta_column, "DeconAllele"]].copy()
            discovery_subset_df.columns = ["index", "discovery beta", "discovery allele"]

            # Save index.
            indices.append("{}\n[N={:,}]".format(discovery_ct, len(ieqtls)))

            ct_replication_data = [np.nan] * len(replication_cell_types)
            ct_annot_data = [""] * len(replication_cell_types)

            if len(ieqtls) > 0:
                ct_replication_data = []
                ct_annot_data = []
                for replication_ct in replication_cell_types:
                    replication_rate = np.nan
                    annot_str = ""

                    # Select the replication eQTLs.
                    replication_subset_df = replication_df.loc[:, ["index", "{} p-value".format(replication_ct), "{} beta".format(replication_ct), "effect_allele"]].copy()
                    replication_subset_df.columns = ["index", "replication p-value", "replication beta", "replication allele"]
                    replication_subset_df.dropna(inplace=True)
                    if replication_subset_df.shape[0] > 0:
                        # Determine overlap.
                        n_overlap = len(set(discovery_subset_df["index"]).intersection(set(replication_subset_df["index"])))

                        # Filter on significant eQTLs.
                        replication_subset_df = replication_subset_df.loc[replication_subset_df["replication p-value"] < 0.05, :]

                        if replication_subset_df.shape[0] > 0:
                            # Merge together.
                            subset_df = discovery_subset_df.merge(replication_subset_df, on=["index"], how="inner")
                            n_replicating = subset_df.shape[0]

                            if subset_df.shape[0] > 0:
                                # Flip direction.
                                subset_df["flip"] = (subset_df["discovery allele"] == subset_df["replication allele"]).replace({False: -1, True: 1})
                                subset_df["replication beta"] = subset_df["replication beta"] * subset_df["flip"]

                                # Calculate concordance.
                                lower_quadrant = subset_df.loc[(subset_df["discovery beta"] < 0) & (subset_df["replication beta"] < 0), :]
                                upper_quadrant = subset_df.loc[(subset_df["discovery beta"] > 0) & (subset_df["replication beta"] > 0), :]
                                concordance = (100 / subset_df.shape[0]) * (lower_quadrant.shape[0] + upper_quadrant.shape[0])

                                # Calculate replication rate.
                                replication_rate = (100 / n_overlap) * n_replicating

                                # Construc the annotation string.
                                annot_str = "Overlap {:,}\nRepl. {} [{:.0f}%]\nConc. {:.0f}%".format(n_overlap, n_replicating, replication_rate, concordance)

                        # Save.
                        ct_replication_data.append(replication_rate)
                        ct_annot_data.append(annot_str)

            # Save.
            replication_data.append(ct_replication_data)
            annot_data.append(ct_annot_data)

        replication_df = pd.DataFrame(replication_data, index=indices, columns=replication_cell_types)
        annot_df = pd.DataFrame(annot_data, index=indices, columns=replication_cell_types)
        print(replication_df)
        print(annot_df)

        # Plot.
        self.plot_heatmap(df=replication_df,
                          annot_df=annot_df)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    def plot_heatmap(self, df, annot_df, xlabel="", ylabel="", appendix=""):
        cmap = sns.diverging_palette(246, 24, as_cmap=True)

        fig, axes = plt.subplots(nrows=2,
                                 ncols=2,
                                 figsize=(1 * df.shape[1] + 10, 1 * df.shape[0] + 10),
                                 gridspec_kw={"width_ratios": [0.2, 0.8],
                                              "height_ratios": [0.8, 0.2]})
        sns.set(color_codes=True)

        annot_df.fillna("", inplace=True)

        row_index = 0
        col_index = 0
        for _ in range(4):
            ax = axes[row_index, col_index]
            if row_index == 0 and col_index == 1:
                sns.heatmap(df, cmap=cmap, vmin=0, vmax=100, center=0,
                            square=True, annot=annot_df, fmt='',
                            cbar=False, annot_kws={"size": 10, "color": "#000000"},
                            ax=ax)

                plt.setp(ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20, rotation=0))
                plt.setp(ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20, rotation=90))

                ax.set_xlabel(xlabel, fontsize=14)
                ax.xaxis.set_label_position('top')

                ax.set_ylabel(ylabel, fontsize=14)
                ax.yaxis.set_label_position('right')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > 1:
                col_index = 0
                row_index += 1

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_corr_heatmap{}.png".format(self.out_filename, appendix)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Discovery:")
        print("    > Data path: {}".format(self.discovery_data_path))
        print("    > Alleles path: {}".format(self.discovery_alleles_path))
        print("    > Name: {}".format(self.discovery_name))
        print("  > Replication:")
        print("    > File path: {}".format(self.replication_path))
        print("    > Name: {}".format(self.replication_name))
        print("  > Output filename: {}".format(self.out_filename))
        print("  > Outpath {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
