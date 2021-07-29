#!/usr/bin/env python3

"""
File:         check_shuffle.py
Created:      2021/07/28
Last Changed: 2021/07/29
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
import glob
import sys
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
./check_shuffle.py -id /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts -if CortexEUR-cis-WithPermutations -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis/data/SampleToDataset.txt.gz
"""

# Metadata
__program__ = "Check Shuffle"
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
        indir = getattr(arguments, 'indir')
        infolder = getattr(arguments, 'infolder')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        self.n_files = getattr(arguments, 'n_files')

        # Set variables.
        if indir is None:
            indir = str(Path(__file__).parent.parent)
        self.indir = os.path.join(indir, "decon_eqtl_with_permutation_fdr", infolder)
        self.outdir = os.path.join(self.indir, "plot")
        for dir in [self.indir, self.outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        self.dataset_cohort_dict = {
            "AMPAD-MAYO-V2": "MAYO",
            "CMC_HBCC_set2": "CMC HBCC",
            "GTEx": "GTEx",
            "AMPAD-ROSMAP-V2": "ROSMAP",
            "BrainGVEX-V2": "Brain GVEx",
            "TargetALS": "Target ALS",
            "AMPAD-MSBB-V2": "MSBB",
            "NABEC-H610": "NABEC",
            "LIBD_1M": "LIBD",
            "ENA": "ENA",
            "LIBD_h650": "LIBD",
            "GVEX": "GVEX",
            "NABEC-H550": "NABEC",
            "CMC_HBCC_set3": "CMC HBCC",
            "UCLA_ASD": "UCLA ASD",
            "CMC": "CMC",
            "CMC_HBCC_set1": "CMC HBCC"
        }

        self.cohort_palette = {
            "MAYO": "#9c9fa0",
            "CMC HBCC": "#0877b4",
            "GTEx": "#0fa67d",
            "ROSMAP": "#6950a1",
            "Brain GVEx": "#48b2e5",
            "Target ALS": "#d5c77a",
            "MSBB": "#5cc5bf",
            "NABEC": "#6d743a",
            "LIBD": "#e49d26",
            "ENA": "#d46727",
            "GVEX": "#000000",
            "UCLA ASD": "#f36d2a",
            "CMC": "#eae453",
            "NA": "#808080"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-id",
                            "--indir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the input path.")
        parser.add_argument("-if",
                            "--infolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the input folder.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=False,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-n",
                            "--n_files",
                            type=int,
                            required=False,
                            default=5,
                            help="The number of files to load. Default: 5.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading color info.")
        col_colors = None
        if self.std_path is not None:
            std_df = self.load_file(self.std_path, header=0, index_col=None)
            col_colors = [self.cohort_palette[self.dataset_cohort_dict[dataset]] for dataset in std_df.iloc[:, 1]]

        ########################################################################

        print("Plotting permutation order.")
        perm_order_m_list = []
        for i, perm_order_inpath in enumerate(glob.glob(os.path.join(self.indir, "perm_orders_*"))):
            for i, perm_pvalues_inpath in enumerate(glob.glob(os.path.join(self.indir, "permutation_pvalues_*"))):
                if i >= (self.n_files - 1):
                    break
            perm_order_m_list.append(self.load_matrix(perm_order_inpath))
        perm_order_m = np.vstack(perm_order_m_list)
        del perm_order_m_list
        print("Shape: {}".format(perm_order_m.shape))

        self.plot_clustermap(m=perm_order_m, col_colors=col_colors,
                             filename="perm_order_m")
        del perm_order_m
        print("")

        ########################################################################

        print("Loading permutation order overlap")
        perm_order_overlap_m_list = []
        for i, perm_order_overlap_inpath in enumerate(glob.glob(os.path.join(self.indir, "perm_order_overlap_*"))):
            for i, perm_pvalues_inpath in enumerate(glob.glob(os.path.join(self.indir, "permutation_pvalues_*"))):
                if i >= (self.n_files - 1):
                    break
            perm_order_overlap_m_list.append(self.load_matrix(perm_order_overlap_inpath))
        perm_order_overlap_m = np.vstack(perm_order_overlap_m_list)
        del perm_order_overlap_m_list
        print("Shape: {}".format(perm_order_overlap_m.shape))

        print("\tplotting distribution")
        self.distplot(x=perm_order_overlap_m.flatten(), xlabel="%overlap", filename="overlap")
        print("")

        ########################################################################

        print("Loading permutation p-value data")
        perm_pvalues_m_list = []
        for i, perm_pvalues_inpath in enumerate(glob.glob(os.path.join(self.indir, "permutation_pvalues_*"))):
            if i >= (self.n_files - 1):
                break
            perm_pvalues_m_list.append(self.load_matrix(perm_pvalues_inpath))
        perm_pvalues_m = np.dstack(perm_pvalues_m_list)
        del perm_pvalues_m_list
        print("Shape: {}".format(perm_pvalues_m.shape))

        print("\tplotting distributions")
        for cov_index in range(perm_pvalues_m.shape[1]):
            self.distplot(x=perm_pvalues_m[:, cov_index, :].flatten(), xlabel="permutation p-value", filename="perm_pvalues_cov{}".format(cov_index))
        print("")

        ########################################################################

        print("Loading real p-value data")
        real_pvalues_df = self.load_file(os.path.join(self.indir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        real_pvalues_m = real_pvalues_df.loc[:, [x for x in real_pvalues_df.columns if x.endswith("_pvalue")]].to_numpy()

        print("\tplotting distribution")
        self.distplot(x=real_pvalues_m.flatten(), xlabel="real p-value", filename="real_pvalues")
        print("")

        ########################################################################

        # Align the shape of the matrices.
        perm_order_overlap_m = np.repeat(perm_order_overlap_m[:, np.newaxis, :], perm_pvalues_m.shape[1], axis=1)
        real_pvalues_m = np.repeat(real_pvalues_m[:, :, np.newaxis], perm_pvalues_m.shape[2], axis=2)

        print("Plotting perm p-value vs %overlap")
        plot_df = pd.DataFrame({"real p-value": real_pvalues_m.flatten(),
                                "perm. p-value": perm_pvalues_m.flatten(),
                                "%overlap": perm_order_overlap_m.flatten()})
        print(plot_df)
        self.regplot(df=plot_df, x="perm. p-value", y="%overlap", group="real p-value", xlabel="perm p-value",
                     ylabel="%identical genotype", filename="pval_vs_overlap")
        print("")

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
    def load_matrix(inpath):
        m = np.load(inpath)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      m.shape))
        return m

    def distplot(self, x, xlabel="", title="", filename="distribution"):
        sns.set_style("ticks")
        fig, ax = plt.subplots(figsize=(12, 12))

        sns.despine(fig=fig, ax=ax)

        sns.kdeplot(x, shade=True, color="#808080", ax=ax, cut=0, zorder=-1)
        ax.axvline(x.mean(), ls='--', color="#808080", zorder=-1)

        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("density",
                      fontsize=14,
                      fontweight='bold')
        ax.set_title(title,
                     fontsize=18,
                     fontweight='bold')

        ax.annotate(
            'N = {:,}'.format(np.size(x)),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            fontsize=15,
            fontweight='bold')

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def regplot(self, df, x, y, group, xlabel="", ylabel="", filename="regplot"):
        offset = 2.2250738585072014e-308

        df.loc[:, [x, group]] = np.log10(df.loc[:, [x, group]] + offset)
        cutoff = np.log10(0.05 + offset)
        df["hue"] = "#808080"
        df.loc[df[x] < cutoff, "hue"] = "#0072B2"
        groups = np.array([-100, -75, -50, -25, -10, -5, -1, 0, None])

        ngroups = np.size(groups)
        ncols = int(np.ceil(np.sqrt(ngroups)))
        nrows = int(np.ceil(ngroups / ncols))

        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=ncols,
                                 sharey="all",
                                 figsize=(12 * ncols, 12 * nrows))
        sns.set(color_codes=True)

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

            if i < ngroups:
                plot_df = df
                title = "All"
                if groups[i] is not None:
                    mask = df[group] <= groups[i]
                    title = "{} <= {}".format(groups[i], group)
                    if i > 0:
                        mask = np.logical_and(df[group] > groups[i - 1], df[group] <= groups[i])
                        title = "{} < {} <= {}".format(groups[i - 1], group, groups[i])
                    plot_df = df.loc[mask, :]

                sns.despine(fig=fig, ax=ax)

                sns.regplot(x=x, y=y, data=plot_df, ci=None,
                            scatter_kws={'facecolors': plot_df["hue"],
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": "#000000"},
                            ax=ax)

                ax.axvline(cutoff, ls='--', color="#000000", zorder=-1)

                ax.annotate(
                    'r = {:.2f}'.format(*stats.spearmanr(plot_df[y], plot_df[x])),
                    xy=(0.03, 0.94),
                    xycoords=ax.transAxes,
                    color="#000000",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'total N = {}'.format(plot_df.shape[0]),
                    xy=(0.03, 0.9),
                    xycoords=ax.transAxes,
                    color="#404040",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'N = {}'.format(plot_df[plot_df["hue"] == "#0072B2"].shape[0]),
                    xy=(0.03, 0.86),
                    xycoords=ax.transAxes,
                    color="#0072B2",
                    fontsize=15,
                    fontweight='bold')
                ax.annotate(
                    'N = {}'.format(plot_df[plot_df["hue"] == "#808080"].shape[0]),
                    xy=(0.03, 0.82),
                    xycoords=ax.transAxes,
                    color="#808080",
                    fontsize=15,
                    fontweight='bold')

                ax.set_xlabel("log10 {}".format(xlabel),
                              fontsize=14,
                              fontweight='bold')
                ax.set_ylabel("{}".format(ylabel),
                              fontsize=14,
                              fontweight='bold')
                ax.set_title(title,
                             fontsize=18,
                             fontweight='bold')
            else:
                ax.set_axis_off()

            col_index += 1
            if col_index > (ncols - 1):
                col_index = 0
                row_index += 1

        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        fig.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def plot_clustermap(self, m, col_colors=None, filename="clustermap"):
        sys.setrecursionlimit(100000)

        sns.set(color_codes=True)
        g = sns.clustermap(m, cmap="Blues",
                           col_colors=col_colors,
                           row_cluster=True, col_cluster=False,
                           yticklabels=False, xticklabels=False,
                           figsize=(12, 9))
        outpath = os.path.join(self.outdir, "{}.png".format(filename))
        plt.tight_layout()
        g.savefig(outpath)
        plt.close()
        print("\tSaved figure: {} ".format(os.path.basename(outpath)))

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Sample-to-dataset path: {}".format(self.std_path))
        print("  > N-files: {}".format(self.n_files))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
