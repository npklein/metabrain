#!/usr/bin/env python3

"""
File:         compare_deconvolution_matrices.py
Created:      2020/09/03
Last Changed: 2021/08/06
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
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Compare Deconvolution Matrices"
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

"""
Syntax:
./compare_deconvolution_matrices.py -d1 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n1 Manuscript -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/cortex_eur_cis/deconvolutionResults.txt.gz -n2 CorrectFTest -log10 -o manuscript_vs_correctFTest

./compare_deconvolution_matrices.py -d1 ../2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -n1 Manuscript -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis/deconvolutionResults.txt.gz -n2 CorrectMDSAndFTest -log10 -o manuscript_vs_correctMDSAndFTest

./compare_deconvolution_matrices.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/cortex_eur_cis/deconvolutionResults.txt.gz -n1 CorrectFTest -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis/deconvolutionResults.txt.gz -n2 CorrectMDSAndFTest -log10 -o correctFTest_vs_correctMDSAndFTest

./compare_deconvolution_matrices.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-WithPermutations-StaticInteractionShuffle/deconvolutionResults.txt.gz -n1 Original -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-HalfNormalizedMAF5/deconvolutionResults.txt.gz -n2 HalfNormalized -log10 -o original_vs_halfnormalized

./compare_deconvolution_matrices.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-NormalisedMAF5/deconvolutionResults.txt.gz -n1 OriginalProfile -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-NormalisedMAF5-CNS7Profile/deconvolutionResults.txt.gz -n2 CNS7Profile -log10 -o originalProfile_vs_CNS7Profile

./compare_deconvolution_matrices.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-NormalisedMAF5-LimitedConfigs/deconvolutionResults.txt.gz -n1 limited_configs -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_with_permutation_fdr/CortexEUR-cis-NormalisedMAF5-ALlConfigs/deconvolutionResults.txt.gz -n2 all_configs -log10 -o limited_vs_all_configurations

./compare_deconvolution_matrices.py -d1 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl/test/deconvolutionResults.txt.gz -n1 original -d2 /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/decon-eqtl_scripts/decon_eqtl_reborn/test/deconvolutionResults.txt.gz -n2 reborn -log10 -o original_vs_reborn
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon1_path = getattr(arguments, 'decon1')
        self.decon1_name = getattr(arguments, 'name1')
        self.decon2_path = getattr(arguments, 'decon2')
        self.decon2_name = getattr(arguments, 'name2')
        self.log10_transform = getattr(arguments, 'log10')
        self.outfile = getattr(arguments, 'outfile')
        self.nrows = None
        self.alpha = 0.05

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'plot')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.accent_colormap = {
            "Intercept": "#000000",
            "Excitatory": "#56B4E9",
            "Inhibitory": "#0072B2",
            "Neuron": "#0072B2",
            "Oligodendrocyte": "#009E73",
            "OPC": "#009E73",
            "EndothelialCell": "#CC79A7",
            "Microglia": "#E69F00",
            "Macrophage": "#E69F00",
            "Astrocyte": "#D55E00",
            "Pericytes": "#808080",
            "Microglia/Macrophage": "#E69F00",
            "Excitatory/Neuron": "#56B4E9",
            "Inhibitory/Neuron": "#0072B2",
            "Excitatory+Inhibitory/Neuron": "#BEBEBE"
        }
        self.point_colormap = {
            "no signif": "#808080",
            "x signif": "#0072B2",
            "y signif": "#D55E00",
            "both signif": "#009E73"
        }

    def create_argument_parser(self):
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-d1",
                            "--decon1",
                            type=str,
                            required=True,
                            help="The path to the first deconvolution matrix")
        parser.add_argument("-n1",
                            "--name1",
                            type=str,
                            required=False,
                            default="x",
                            help="The name for the first deconvolution matrix")
        parser.add_argument("-d2",
                            "--decon2",
                            type=str,
                            required=True,
                            help="The path to the second deconvolution matrix")
        parser.add_argument("-n2",
                            "--name2",
                            type=str,
                            required=False,
                            default="y",
                            help="The name for the first deconvolution matrix")
        parser.add_argument("-log10",
                            action='store_true',
                            help="Log10 transform the profile values."
                                 " Default: False.")
        parser.add_argument("-o",
                            "--outfile",
                            type=str,
                            required=False,
                            default="deconvolution",
                            help="The name of the plot file.")

        return parser.parse_args()

    def start(self):

        print("Loading data")
        decon1_df = self.load_file(self.decon1_path)
        decon2_df = self.load_file(self.decon2_path)

        print("Splitting categories")
        pval_df1, fdr_df1, intercept_beta_df1, ct_beta_df1, inter_beta_df1 = self.split_decon(df=decon1_df)
        pval_df2, fdr_df2, intercept_beta_df2, ct_beta_df2, inter_beta_df2 = self.split_decon(df=decon2_df)

        # inter_beta_df1 = np.abs(inter_beta_df1)
        # inter_beta_df2 = np.abs(inter_beta_df2)

        print("Plotting")
        for decon1_df, decon2_df, appendix, log10_transform, signif_lines in \
                ([intercept_beta_df1, intercept_beta_df2, "InterceptBeta", False, False],
                 [ct_beta_df1, ct_beta_df2, "CellTypeBeta", False, False],
                 [inter_beta_df1, inter_beta_df2, "IteractionBeta", False, False],
                 [pval_df1, pval_df2, "pvalue", self.log10_transform, True],
                 [fdr_df1, fdr_df2, "FDR", self.log10_transform, True]):
            if decon1_df is None and decon2_df is None:
                continue

            print("Plotting {}".format(appendix))
            self.plot_violin_comparison(decon1_df=decon1_df,
                                        decon2_df=decon2_df,
                                        filename=self.outfile + "_" + appendix)

            if decon1_df is None or decon2_df is None:
                continue

            # Filter on overlapping rows and columns.
            print("\tFiltering matrices.")
            row_overlap = set(decon1_df.index).intersection(set(decon2_df.index))
            decon1_df = decon1_df.loc[row_overlap, :]
            decon2_df = decon2_df.loc[row_overlap, :]

            row_overlap = set(decon1_df.index).intersection(set(decon2_df.index))
            col_overlap = set(decon1_df.columns).intersection(set(decon2_df.columns))

            decon1_df = decon1_df.loc[row_overlap, col_overlap]
            decon2_df = decon2_df.loc[row_overlap, col_overlap]

            print("\t  Deconvolution matrix 1: {}".format(decon1_df.shape))
            print("\t  Deconvolution matrix 2: {}".format(decon2_df.shape))

            if decon1_df.shape != decon2_df.shape:
                print("\t  Shape's are not identical.")
                exit()

            decon1_signif_df = None
            decon2_signif_df = None
            if "Beta" in appendix:
                decon1_signif_df = fdr_df1
                decon2_signif_df = fdr_df2

            print("\tPlotting.")
            self.plot_regression_comparison(decon1_df=decon1_df,
                                            decon2_df=decon2_df,
                                            decon1_signif_df=decon1_signif_df,
                                            decon2_signif_df=decon2_signif_df,
                                            log10_transform=log10_transform,
                                            signif_lines=signif_lines,
                                            filename=self.outfile + "_" + appendix)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def split_decon(df):
        # Get the p-values.
        pval_df = df.loc[:, [x for x in df.columns if x.endswith("_pvalue")]].copy()
        pval_df.columns = [x.split("_")[1] for x in pval_df]

        # Calculate the FDR.
        fdr_df = pd.DataFrame(np.nan, columns=pval_df.columns, index=pval_df.index)
        for index, row in pval_df.T.iterrows():
            fdr_df[index] = multitest.multipletests(row, method='fdr_bh')[1]

        # Intercept beta.
        intercept_beta_df = None
        if "Beta1_Intercept" in df.columns:
            intercept_beta_df = df.loc[:, ["Beta1_Intercept"]]
            intercept_beta_df.columns = ["Intercept"]

        # Get the cell type beta's.
        ct_beta_df = df.loc[:, [x for x in df.columns if "Beta" in x and not ":GT" in x and not x in ["Beta1_Intercept", "Beta2_Genotype"]]]
        ct_beta_df.columns = [x.split("_")[2] for x in ct_beta_df.columns]

        # Get the interaction beta's.
        interaction_beta_df = df.loc[:, [x for x in df.columns if "Beta" in x and ":GT" in x]]
        interaction_beta_df.columns = [x.split("_")[2].split(":")[0] for x in interaction_beta_df.columns]

        return pval_df, fdr_df, intercept_beta_df, ct_beta_df, interaction_beta_df

    def print_n_signif(self, df, type):
        print("\nN-interaction ({} < {}):".format(type, self.alpha))
        n_hits_s = (df < self.alpha).sum(axis=0)
        n_hits_total = n_hits_s.sum()
        cov_length = np.max([len(x) for x in df.columns])
        hits_length = np.max([len(str(x)) for x in n_hits_s] + [len(str(n_hits_total))])
        for cell_type, n_hits in n_hits_s.iteritems():
            print("\t{:{}s}  {:{}d}".format(cell_type, cov_length, n_hits, hits_length))
        print("\t{}".format("".join(["-"] * cov_length)))
        print("\t{:{}s}  {:{}d}".format("total", cov_length, n_hits_total, hits_length))

        print("", flush=True)

    def plot_violin_comparison(self, decon1_df, decon2_df, filename):
        df1 = None
        if decon1_df is not None:
            df1 = decon1_df.copy()
            df1.columns = [x.replace(":GT", "") for x in df1.columns]
            df1 = df1.melt()
            df1["name"] = self.decon1_name

        df2 = None
        if decon2_df is not None:
            df2 = decon2_df.copy()
            df2.columns = [x.replace(":GT", "") for x in df2.columns]
            df2 = df2.melt()
            df2["name"] = self.decon2_name

        df = pd.concat([df1, df2], axis=0)

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)
        sns.violinplot(x="name", y="value", hue="variable", legend=False,
                       data=df, palette=self.accent_colormap, ax=ax)
        ax.set_xlabel("",
                      fontsize=14,
                      fontweight='bold')
        ax.set_ylabel("value",
                      fontsize=14,
                      fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_violin_comparison.png".format(filename)))
        plt.close()

    def plot_regression_comparison(self, decon1_df, decon2_df, decon1_signif_df=None,
                                   decon2_signif_df=None, log10_transform=False,
                                   signif_lines=False, filename="plot"):

        # celltypes = (("Microglia", "Macrophage"),
        #              ("EndothelialCell", "EndothelialCell"),
        #              ("Oligodendrocyte", "Oligodendrocyte"),
        #              ("Excitatory", "Neuron"),
        #              ("Inhibitory", "Neuron"),
        #              ("Astrocyte", "Astrocyte"))

        # Calculate number of rows / columns.
        celltypes = list(decon1_df.columns)
        celltypes.sort()
        ncols = 3
        nrows = math.ceil(len(celltypes) / 3)

        sns.set(rc={'figure.figsize': (12 * ncols, 9 * nrows)})
        sns.set_style("ticks")
        fig = plt.figure()
        grid = fig.add_gridspec(ncols=ncols,
                                nrows=nrows)

        row_index = 0
        col_index = 0

        # for ct2, ct1 in celltypes:
        for celltype in celltypes:
            x_label = "{}_{}".format(self.decon1_name, celltype)
            y_label = "{}_{}".format(self.decon2_name, celltype)
            df = pd.DataFrame({x_label: decon1_df.loc[:, celltype],
                               y_label: decon2_df.loc[:, celltype]})

            x_signif_col = x_label
            y_signif_col = y_label
            if decon1_signif_df is not None and decon2_signif_df is not None:
                df["x_signif"] = decon1_signif_df.loc[:, celltype]
                df["y_signif"] = decon2_signif_df.loc[:, celltype]
                x_signif_col = "x_signif"
                y_signif_col = "y_signif"

            # Add color
            df["hue"] = self.point_colormap["no signif"]
            df.loc[(df[x_signif_col] < 0.05) & (df[y_signif_col] >= 0.05), "hue"] = self.point_colormap["x signif"]
            df.loc[(df[x_signif_col] >= 0.05) & (df[y_signif_col] < 0.05), "hue"] = self.point_colormap["y signif"]
            df.loc[(df[x_signif_col] < 0.05) & (df[y_signif_col] < 0.05), "hue"] = self.point_colormap["both signif"]

            counts = df["hue"].value_counts()
            for value in self.point_colormap.values():
                if value not in counts.index:
                    counts[value] = 0

            coef, _ = stats.spearmanr(df[x_label], df[y_label])
            coef_str = "{:.2f}".format(coef)

            accent_color = self.accent_colormap[celltype]

            # Convert to log10 scale.
            if log10_transform:
                df[x_label] = np.log10(df[x_label] + 1)
                df[y_label] = np.log10(df[y_label] + 1)

            # Creat the subplot.
            ax = fig.add_subplot(grid[row_index, col_index])
            sns.despine(fig=fig, ax=ax)

            # Plot.
            g = sns.regplot(x=x_label,
                            y=y_label,
                            data=df,
                            scatter_kws={'facecolors': df["hue"],
                                         'linewidth': 0,
                                         'alpha': 0.5},
                            line_kws={"color": accent_color},
                            ax=ax
                            )

            # sns.jointplot(data=df,
            #               x=x_label,
            #               y=y_label,
            #               kind="reg",
            #               scatter_kws={'facecolors': df["hue"],
            #                            'linewidth': 0,
            #                            'alpha': 0.5},
            #               line_kws={"color": accent_color}
            #               )

            # Add the text.
            ax.annotate(
                'r = {}'.format(coef_str),
                xy=(0.03, 0.94),
                xycoords=ax.transAxes,
                color="#404040",
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'total N = {:,}'.format(df.shape[0]),
                xy=(0.03, 0.9),
                xycoords=ax.transAxes,
                color="#404040",
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[self.point_colormap["both signif"]]),
                xy=(0.03, 0.86),
                xycoords=ax.transAxes,
                color=self.point_colormap["both signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[self.point_colormap["x signif"]]),
                xy=(0.03, 0.82),
                xycoords=ax.transAxes,
                color=self.point_colormap["x signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[self.point_colormap["y signif"]]),
                xy=(0.03, 0.78),
                xycoords=ax.transAxes,
                color=self.point_colormap["y signif"],
                fontsize=18,
                fontweight='bold')
            ax.annotate(
                'N = {:,}'.format(counts[self.point_colormap["no signif"]]),
                xy=(0.03, 0.74),
                xycoords=ax.transAxes,
                color=self.point_colormap["no signif"],
                fontsize=18,
                fontweight='bold')

            ax.text(0.5, 1.05,
                    celltype,
                    fontsize=26, weight='bold', ha='center', va='bottom',
                    color=accent_color,
                    transform=ax.transAxes)

            ax.set_xlabel(" ".join(x_label.split("_")),
                          fontsize=18,
                          fontweight='bold')

            ax.set_ylabel(" ".join(y_label.split("_")),
                          fontsize=18,
                          fontweight='bold')

            if signif_lines:
                signif_line = 0.05
                if log10_transform:
                    signif_line = np.log10(1.05)
                ax.axhline(signif_line, ls='--', color="#000000", zorder=-1)
                ax.axvline(signif_line, ls='--', color="#000000", zorder=-1)

            # min_value = np.min((df[[x_label, y_label]].min(axis=0).min() * 1.1, -0.1))
            # max_value = np.max((df[[x_label, y_label]].max(axis=0).max() * 1.1, 1.1))
            #
            # print(min_value, max_value)
            #
            # ax.set_xlim(min_value, max_value)
            # ax.set_ylim(min_value, max_value)
            # ax.plot([min_value, max_value], [min_value, max_value], ls="--", c="#000000")

            # Increment indices.
            col_index += 1
            if col_index == ncols:
                col_index = 0
                row_index += 1

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_regression_comparison.png".format(filename)))
        plt.close()


if __name__ == '__main__':
    m = main()
    m.start()
