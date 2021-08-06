#!/usr/bin/env python3

"""
File:         visualise_permutation_interactions.py
Created:      2021/08/05
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
import warnings
import argparse
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.

# Metadata
__program__ = "Visualise Permutation Interactions"
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
./visualise_permutation_interactions.py -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/select_and_reorder_matrix/CortexEUR-cis/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/perform_deconvolution/deconvolution_table.txt.gz -std /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_cohort_matrix/sample_to_dataset.txt.gz -df CortexEUR-cis-WithPermutations-StaticGenoShuffle
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.std_path = getattr(arguments, 'sample_to_dataset')
        decon_folder = getattr(arguments, 'decon_folder')
        self.n = getattr(arguments, 'n_plots')

        # Set variables.
        self.decondir = os.path.join(str(Path(__file__).parent.parent), "decon_eqtl_with_permutation_fdr", decon_folder)
        self.outdir = os.path.join(self.decondir, "plot")
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.colormap = {
            "minor": "#E69F00",
            "center": "#0072B2",
            "major": "#D55E00"
        }

        # Create color map.
        self.palette = {0.0: "#D55E00",
                        1.0: "#0072B2",
                        2.0: "#E69F00"}

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
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix.")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix.")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts matrix.")
        parser.add_argument("-std",
                            "--sample_to_dataset",
                            type=str,
                            required=True,
                            help="The path to the sample-to-dataset matrix.")
        parser.add_argument("-df",
                            "--decon_folder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the deconvolution folder.")
        parser.add_argument("-n",
                            "--n_plots",
                            type=int,
                            required=False,
                            default=5,
                            help="The number of plots to create. Default: 5.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading permutation order.")
        perm_order_inpath = glob.glob(os.path.join(self.decondir, "perm_orders_*"))
        perm_order_inpath.sort()
        perm_order_m_list = []
        for perm_order_inpath in perm_order_inpath:
            perm_order_m_list.append(self.load_matrix(perm_order_inpath))

        perm_order_m = np.vstack(perm_order_m_list)
        del perm_order_m_list
        print("\tShape: {}".format(perm_order_m.shape))

        print("Loading permutation p-value data")
        perm_pvalues_inpaths = glob.glob(os.path.join(self.decondir, "permutation_pvalues_*"))
        perm_pvalues_inpaths.sort()
        perm_pvalues_m_list = []
        for perm_pvalues_inpath in perm_pvalues_inpaths:
            perm_pvalues_m_list.append(self.load_matrix(perm_pvalues_inpath))
        perm_pvalues_m = np.dstack(perm_pvalues_m_list)
        del perm_pvalues_m_list
        print("\tShape: {}".format(perm_pvalues_m.shape))

        print("Loading genotype data")
        geno_m = self.load_file(self.geno_path, header=0, index_col=0).to_numpy(dtype=np.float64)

        print("Calculate MAF")
        maf_a = np.empty(geno_m.shape[0], dtype=np.float64)
        for geno_index in range(geno_m.shape[0]):
            unique, counts = np.unique(geno_m[geno_index, :], return_counts=True)
            group_sizes = dict(zip(unique, counts))
            group_sizes_a = np.array([group_sizes.get(0, 0),
                                      group_sizes.get(1, 0),
                                      group_sizes.get(2, 0)],
                                     dtype=np.float64)

            allele1 = group_sizes_a[0] * 2 + group_sizes_a[1]
            allele2 = group_sizes_a[2] * 2 + group_sizes_a[1]
            maf_a[geno_index] = min(allele1, allele2) / (allele1 + allele2)
        maf_mask = maf_a > 0.14

        print("Finding the lowest N values.")
        maf_filtered_perm_pvalues_m = perm_pvalues_m[maf_mask, :]
        smallest_indices = []
        cutoff = 0
        nrows = None
        for i in range(self.n):
            lowest_value = np.min(maf_filtered_perm_pvalues_m[maf_filtered_perm_pvalues_m > cutoff])
            index = np.where(perm_pvalues_m == lowest_value)

            row_index = index[0][0]
            col_index = index[1][0]
            perm_index = index[2][0]

            smallest_indices.append((row_index, col_index, perm_index))
            cutoff = lowest_value

            if nrows is None or row_index > nrows:
                nrows = row_index

        print(smallest_indices)

        print("Loading real p-value data")
        real_pvalues_df = self.load_file(os.path.join(self.decondir, "deconvolutionResults.txt.gz"), header=0, index_col=0)
        real_pvalues_m = real_pvalues_df.loc[:, [x for x in real_pvalues_df.columns if x.endswith("_pvalue")]].to_numpy()

        print("Loading other data")
        geno_m = self.load_file(self.geno_path, header=0, index_col=0, nrows=nrows + 1).to_numpy(dtype=np.float64)
        expr_m = self.load_file(self.expr_path, header=0, index_col=0, nrows=nrows + 1).to_numpy(dtype=np.float64)
        cc_m = self.load_file(self.cc_path, header=0, index_col=0).to_numpy(dtype=np.float64)
        std_m = self.load_file(self.std_path, header=0, index_col=None).to_numpy(dtype=object)

        print("Filling genotype nan with dataset mean")
        geno_m[geno_m == -1] = np.nan

        dataset_counts = list(zip(*np.unique(std_m[:, 1], return_counts=True)))
        dataset_counts.sort(key=lambda x: -x[1])
        datasets = [x[0] for x in dataset_counts]

        geno_nan_mask = np.isnan(geno_m)
        nanfilled_geno_m = np.copy(geno_m)
        if np.sum(geno_nan_mask) > 0:
            # Calculate the dataset genotype means per eQTL and fill nan values
            # with the dataset mean.
            geno_dataset_mean_m = self.calculate_geno_mean_per_dataset(geno_m=geno_m,
                                                                       datasets=datasets,
                                                                       std_m=std_m)
            nanfilled_geno_m[geno_nan_mask] = geno_dataset_mean_m[geno_nan_mask]

        print("Plotting")
        for row_index, cov_index, perm_index in smallest_indices:
            genotype = geno_m[row_index, :]
            expression = expr_m[row_index, :]
            cell_count = cc_m[:, cov_index]
            perm_order = perm_order_m[perm_index, :]
            nanfilled_genotype = nanfilled_geno_m[row_index, :]
            shuffled_genotype = nanfilled_genotype[perm_order]

            plot_df = pd.DataFrame({"genotype": genotype,
                                    "genotype group": np.rint(genotype),
                                    "expression": expression,
                                    "cell count": cell_count,
                                    "shuffled genotype": shuffled_genotype,
                                    "shuffled genotype group": np.rint(shuffled_genotype)})

            real_p_value = real_pvalues_m[row_index, cov_index]
            perm_p_value = perm_pvalues_m[row_index, cov_index, perm_index]

            # Initialize plot.
            sns.set(rc={'figure.figsize': (24, 18)})
            sns.set_style("ticks")
            fig, axes = plt.subplots(ncols=2, nrows=2)

            # Plot the original genotype.
            self.plot_eqtl(df=plot_df.loc[plot_df['genotype'].notnull(), :],
                           x="genotype",
                           y="expression",
                           x_group="genotype group",
                           palette=self.palette,
                           ax=axes[0, 0],
                           title="real genotype")
            self.plot_inter_eqtl(df=plot_df.loc[plot_df['genotype'].notnull(), :],
                                 x="cell count",
                                 y="expression",
                                 group="genotype group",
                                 palette=self.palette,
                                 ax=axes[1, 0],
                                 annotate=[("Decon-eQTL p-value", real_p_value, ".2e")])

            self.plot_eqtl(df=plot_df.loc[plot_df['shuffled genotype'].notnull(), :],
                           x="shuffled genotype",
                           y="expression",
                           x_group="shuffled genotype group",
                           palette=self.palette,
                           ax=axes[0, 1],
                           title="permuted genotype")
            self.plot_inter_eqtl(df=plot_df.loc[plot_df['shuffled genotype'].notnull(), :],
                                 x="cell count",
                                 y="expression",
                                 group="shuffled genotype group",
                                 palette=self.palette,
                                 ax=axes[1, 1],
                                 annotate=[("Decon-eQTL p-value", perm_p_value, ".2e")])

            outpath = os.path.join(self.outdir, "permutation_interaction_{}_{}_{}.png".format(row_index, cov_index, perm_index))
            fig.savefig(outpath)
            plt.close()
            print("\tSaved: {}".format(outpath))

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=None, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def load_matrix(inpath):
        m = np.load(inpath)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      m.shape))
        return m

    @staticmethod
    def calculate_geno_mean_per_dataset(geno_m, datasets, std_m):
        """
        Method for calculating the mean genotype per dataset per eQTL. Missing
        values are not included. Returns a matrix with n-eQTLs x n-samples
        with the mean genotype per dataset in the cells.
        """
        geno_dataset_mean_m = np.empty_like(geno_m)

        # Ignore RuntimeWarning: Mean of empty slice
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            # Fill in the missing values with the dataset mean.
            for dataset_index, dataset in enumerate(datasets):
                dataset_mask = std_m[:, 1] == dataset
                geno_dataset_m = geno_m[:, dataset_mask]
                dataset_mean_a = np.nanmean(geno_dataset_m, axis=1)
                geno_dataset_mean_m[:, dataset_mask] = np.tile(dataset_mean_a[:, np.newaxis], np.sum(dataset_mask))

        return geno_dataset_mean_m

    @staticmethod
    def plot_eqtl(df, x, y, x_group, palette, ax, title="", xlabel="SNP", ylabel="gene",
                  annotate=None):
        # Calculate the correlation.
        coef, _ = stats.pearsonr(df[x], df[y])

        unique, counts = np.unique(df[x_group].to_numpy(), return_counts=True)
        group_sizes = dict(zip(unique, counts))
        group_sizes_a = np.array([group_sizes.get(0, 0),
                                  group_sizes.get(1, 0),
                                  group_sizes.get(2, 0)],
                                 dtype=np.float64)

        allele1 = group_sizes_a[0] * 2 + group_sizes_a[1]
        allele2 = group_sizes_a[2] * 2 + group_sizes_a[1]
        maf = min(allele1, allele2) / (allele1 + allele2)

        # Plot the scatter / box plot.
        sns.regplot(x=x, y=y, data=df,
                    scatter=False,
                    line_kws={"color": "#000000"},
                    ax=ax
                    )
        sns.boxplot(x=x_group, y=y, data=df,
                    palette=palette,
                    showfliers=False,
                    zorder=1,
                    ax=ax)

        ax.annotate(
            'N = {:,}'.format(df.shape[0]),
            xy=(0.03, 0.94),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=12,
            fontweight='bold')
        ax.annotate(
            'r = {:.2f}'.format(coef),
            xy=(0.03, 0.90),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=12,
            fontweight='bold')
        ax.annotate(
            'MAF = {:.2f}'.format(maf),
            xy=(0.03, 0.86),
            xycoords=ax.transAxes,
            color="#000000",
            alpha=1,
            fontsize=12,
            fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.86 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=12,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    @staticmethod
    def plot_inter_eqtl(df, x, y, group, palette, ax, title="SNP",
                        xlabel="gene", ylabel="", annotate=None):

        for i, genotype in enumerate([0.0, 1.0, 2.0]):
            subset = df.loc[df[group] == genotype, :].copy()
            color = palette[genotype]
            coef = np.nan
            if len(subset.index) > 1:
                # Calculate the correlation.
                coef, _ = stats.pearsonr(df[x], df[y])

                # Plot the scatter / box plot.
                sns.regplot(x=x, y=y, data=subset,
                            scatter_kws={'facecolors': color,
                                         'linewidth': 0,
                                         'alpha': 0.75},
                            line_kws={"color": color, "alpha": 0.75},
                            ax=ax
                            )

            ax.annotate(
                '{}: r = {:.2f} [N = {:,}]'.format(genotype, coef, subset.shape[0]),
                xy=(0.03, 0.94 - (0.04 * i)),
                xycoords=ax.transAxes,
                color=color,
                alpha=0.75,
                fontsize=12,
                fontweight='bold')

        if annotate is not None:
            for i, (label, value, rounding) in enumerate(annotate):
                ax.annotate(
                    '{} = {:{}}'.format(label, value, rounding),
                    xy=(0.03, 0.82 - (0.04 * i)),
                    xycoords=ax.transAxes,
                    color="#000000",
                    alpha=0.75,
                    fontsize=12,
                    fontweight='bold')

        ax.set_title(title,
                     fontsize=22,
                     weight='bold')
        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Decon-eQTL directory: {}".format(self.decondir))
        print("  > N-plots: {}".format(self.n))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
