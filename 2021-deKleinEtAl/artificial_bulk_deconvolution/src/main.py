"""
File:         main.py
Created:      2020/11/09
Last Changed: 2020/11/12
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
import random
import glob
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy import stats
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Local application imports.
from .utilities import prepare_output_dir, load_dataframe


class Main:
    def __init__(self, input_path, input_suffix, cell_counts_path,
                 ref_profile_path, gene_info_path, combine, n_samples):
        self.input_path = input_path
        self.input_suffix = input_suffix
        self.cell_counts_path = cell_counts_path
        self.ref_profile_path = ref_profile_path
        self.gene_info_path = gene_info_path
        self.combine = combine
        self.extension = "png"

        if self.input_suffix is None:
            self.input_suffix = ""

        # Define the working directory.
        work_dir = str(Path(__file__).parent.parent)

        # Prepare an output directory.
        self.outdir = os.path.join(work_dir, 'output')
        prepare_output_dir(self.outdir)

        self.cell_type_list = [("END", "END", "Endothelial Cell", "#CC79A7"),
                               ("MIC", "MAC", "Microglia VS Macrophage", "#E69F00"),
                               ("AST", "AST", "Astrocyte", "#D55E00"),
                               ("IN", "NEU", "In. Neuron VS Neuron", "#0072B2"),
                               ("OLI", "OLI", "Oligodendrocyte", "#009E73"),
                               ("EX", "NEU", "Ex. Neuron VS Neuron", "#0072B2")]

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        ref_profile_df = self.load_reference_profile()
        ref_profile_df = self.filter(ref_profile_df, 5)
        ref_profile_df = self.perform_log2_transform(ref_profile_df)
        if ref_profile_df.values.min() < 0:
            ref_profile_df = self.perform_shift(ref_profile_df)

        print("")
        print("### Step2 ###")
        samples, matrices = self.load_ct_matrices(genes=ref_profile_df.index)

        print("")
        print("### Step3 ###")
        cc_df, cell_types = self.load_cell_counts()
        print(cc_df)

        if len(set(list(matrices.keys())).symmetric_difference(set(cell_types))) != 0:
            print("Missing cell type info")
            exit()

        print("")
        print("### Step4 ###")
        sample_matrices = self.merge_samples(samples=samples,
                                             cell_types=cell_types,
                                             matrices=matrices)

        print("")
        print("### Step5 ###")
        real_data = []
        predict_data = []
        indices = []
        for sample in samples:
            if sample not in cc_df:
                continue
            real_weights = cc_df.loc[:, sample]

            _, bulk_expression = self.create_artifical_bulk_data(df=sample_matrices[sample],
                                                                 weights=real_weights)

            # real_weights, bulk_expression = self.create_artifical_bulk_data(df=sample_matrices[sample])

            decon_weights = self.deconvolute(ref_profile_df, bulk_expression)

            # weights_df = pd.DataFrame({"real": real_weights,
            #                            "predict": decon_weights})
            # weights_df.index = ref_profile_df.columns
            # weights_df["diff"] = (weights_df["real"] - weights_df["predict"]).abs()
            # print(weights_df)

            real_data.append(real_weights)
            predict_data.append(decon_weights)
            indices.append(sample)

        real_weights_df = pd.DataFrame(real_data, columns=cell_types, index=indices)
        print(real_weights_df)
        decon_weights_df = pd.DataFrame(predict_data, columns=ref_profile_df.columns, index=indices)
        print(decon_weights_df)

        print("")
        print("### Step6 ###")
        print("Visualizing differences")

        self.plot_spider(real_weights_df, decon_weights_df, cell_types)
        self.plot_regression(real_weights_df, decon_weights_df, cell_types)

        real_weights_df_m = real_weights_df.melt()
        real_weights_df_m["hue"] = "real"
        decon_weights_df_m = decon_weights_df.melt()
        decon_weights_df_m["hue"] = "predict"

        dfm = pd.concat([real_weights_df_m, decon_weights_df_m], axis=0)
        self.plot_distributions(dfm)

    def load_reference_profile(self):
        print("Load reference profile")
        ref_profile_df = load_dataframe(self.ref_profile_path,
                                        sep="\t",
                                        header=0,
                                        index_col=0)

        gene_trans_df = load_dataframe(self.gene_info_path,
                                       sep="\t",
                                       header=0,
                                       index_col=0)[["Symbol", "ArrayAddress"]]
        gene_trans_df.dropna(axis=0, inplace=True)
        gene_trans_dict = dict(zip(gene_trans_df.iloc[:, 0], gene_trans_df.iloc[:, 1]))

        mask = []
        names = []
        for val in ref_profile_df.index:
            if val in gene_trans_dict:
                mask.append(True)
                names.append(gene_trans_dict[val])
            else:
                mask.append(False)

        ref_profile_df = ref_profile_df.loc[mask, :]
        ref_profile_df.index = names

        return ref_profile_df

    @staticmethod
    def filter(df, cutoff):
        print("\tFiltering uninformative genes from the signature matrix")
        return df.loc[(df.std(axis=1) != 0) & (df.max(axis=1) > cutoff), :]

    @staticmethod
    def perform_log2_transform(df):
        print("\tPerforming log2 transformation")
        df = df + abs(df.values.min()) + 1
        return np.log2(df)

    @staticmethod
    def perform_shift(df):
        return df + abs(df.values.min())

    def load_ct_matrices(self, genes):
        print("Loading cell type matrices")
        columns = None
        matrices = {}

        for inpath in glob.glob(os.path.join(self.input_path, "*" + self.input_suffix)):
            cell_type = os.path.basename(inpath).replace(self.input_suffix, "")

            df = load_dataframe(inpath,
                                sep="\t",
                                header=0,
                                index_col=0)
            if df.values.min() < 0:
                df = self.perform_shift(df)
            df = self.perform_log2_transform(df)

            if columns is None:
                columns = df.columns
            else:
                if not columns.equals(df.columns):
                    print("Columns of cell type expression are not identical")
                    exit()

            matrices[cell_type] = df

        full_matrices = {}
        for key, df in matrices.items():
            missing = []
            for gene in genes:
                if gene not in df.index:
                    missing.append(gene)
            if len(missing) > 0:
                missing_df = pd.DataFrame(0, index=missing, columns=columns)
                full_df = df.T.merge(missing_df.T, left_index=True, right_index=True).T
                full_df = full_df.loc[genes, :]
            else:
                full_df = df
            full_matrices[key] = full_df

        del matrices

        return columns, full_matrices

    def load_cell_counts(self):
        df = load_dataframe(self.cell_counts_path,
                            sep="\t",
                            header=0,
                            index_col=0)

        ratios_df = df / df.sum(axis=0)
        order_df = df.mean(axis=1)
        order_df.sort_values(inplace=True, ascending=True)
        order = list(order_df.index)

        return ratios_df.loc[order, :], order

    @staticmethod
    def merge_samples(samples, cell_types, matrices):
        sample_matrices = {}

        for i, sample in enumerate(samples):
            print("\t[{}] Combining cell type expression for sample {}".format(i, sample))
            sample_df = None

            for ct in cell_types:
                sample_ct_values = matrices[ct][[sample]]
                sample_ct_values.columns = [ct]
                if sample_df is None:
                    sample_df = sample_ct_values
                else:
                    sample_df = sample_df.merge(sample_ct_values,
                                                left_index=True,
                                                right_index=True)

            sample_matrices[sample] = sample_df

        return sample_matrices

    @staticmethod
    def create_artifical_bulk_data(df, weights=None):
        if weights is None:
            weights = np.array([random.uniform(0, 1) for _ in range(df.shape[1])])
            weights = weights / np.sum(weights)

        return weights, pd.Series(np.dot(df, weights), index=df.index)

    @staticmethod
    def deconvolute(profile_df, bulk_array):
        weights = nnls(profile_df, bulk_array)[0]
        #return weights / np.sum(weights)
        return weights

    def plot_spider(self, real_df, decon_df, cell_types):
        my_dpi = 96
        fig = plt.figure(figsize=(18, 6), dpi=my_dpi)

        # Loop to plot
        for i, (sn_ct, decon_ct, title, color) in enumerate(self.cell_type_list):
            real = real_df[[sn_ct]]
            real.columns = ["real"]
            decon = decon_df[[decon_ct]]
            decon.columns = ["predict"]
            df = real.merge(decon, left_index=True,right_index=True)
            df.sort_values(by="real", ascending=False, inplace=True)

            N = df.shape[0] - 1
            angles = [n / float(N) * 2 * np.pi for n in range(N)]
            angles += angles[:1]

            ax = plt.subplot(1, len(cell_types), i+1, polar=True)

            ax.set_theta_offset(np.pi / 2)
            ax.set_theta_direction(-1)

            plt.xticks(angles[:-1], [""] * len(angles[:-1]), color='white', size=8)

            ax.set_rlabel_position(0)
            ticks = []
            for x in np.arange(0, 1.25, 0.25):
                if (df.values.max() + 0.25) > x:
                    ticks.append(x)
            plt.yticks(ticks, ticks, color="grey", size=7)
            plt.ylim(0, max(ticks))

            ax.plot(angles, df["real"], color="#808080", linewidth=2, linestyle='solid')
            ax.plot(angles, df["predict"], color=color, linewidth=2, linestyle='solid')

            # Add a title
            plt.title(title, size=11, y=1.1)

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "weights_radarplot.{}".format(self.extension)))
        plt.close()

    def plot_regression(self, real_df, decon_df, cell_types):
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=1,
                                 ncols=len(self.cell_type_list),
                                 figsize=(8*len(self.cell_type_list), 6))

        # Loop to plot
        for i, (sn_ct, decon_ct, title, color) in enumerate(
                self.cell_type_list):
            real = real_df[[sn_ct]]
            real.columns = ["real"]
            decon = decon_df[[decon_ct]]
            decon.columns = ["predict"]
            df = real.merge(decon, left_index=True, right_index=True)

            include_ylabel = False
            if i == 0:
                include_ylabel = True

            self.regplot(df=df,
                         fig=fig,
                         ax=axes[i],
                         x="real",
                         y="predict",
                         xlabel="real cc%",
                         ylabel="predicted cc%",
                         title=title,
                         color=color,
                         include_ylabel=include_ylabel)

        fig.savefig(os.path.join(self.outdir, "weights_correlation.{}".format(self.extension)))
        plt.close()

    def regplot(self, df, fig, ax, x="x", y="y", facecolors=None,
             xlabel="", ylabel="", title="", color="#000000",
             include_ylabel=True):
        sns.despine(fig=fig, ax=ax)

        if not include_ylabel:
            ylabel = ""

        if facecolors is None:
            facecolors = "#808080"
        else:
            facecolors = df[facecolors]

        n = df.shape[0]
        coef = np.nan

        if n > 0:
            coef, p = stats.spearmanr(df[x], df[y])
            # coef, p = stats.pearsonr(df[x], df[y])

            sns.regplot(x=x, y=y, data=df,
                        scatter_kws={'facecolors': facecolors,
                                     'edgecolors': "#808080"},
                        line_kws={"color": color},
                        ax=ax
                        )

        ax.text(0.5, 1.1, title,
                fontsize=18, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02, "N = {}".format(n),
                fontsize=14, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel(ylabel,
                      fontsize=14,
                      fontweight='bold')
        ax.set_xlabel(xlabel,
                      fontsize=14,
                      fontweight='bold')

        ax.legend(handles=[mpatches.Patch(color=color, label="r = {:.2f}".format(coef))], loc=4)

    def plot_distributions(self, df, name=""):
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=1,
                                 ncols=len(self.cell_type_list),
                                 figsize=(8*len(self.cell_type_list), 6))

        for i, (sn_ct, decon_ct, title, color) in enumerate(self.cell_type_list):
            ax = axes[i]
            sns.despine(fig=fig, ax=ax)

            real_values = df.loc[(df["variable"] == sn_ct) & (df["hue"] == "real"), "value"]
            real_values.name = "real"
            predict_values = df.loc[(df["variable"] == decon_ct) & (df["hue"] == "predict"), "value"]
            predict_values.name = "predict"

            if max(real_values) > 0:
                try:
                    sns.kdeplot(real_values, shade=True, color="#808080", ax=ax,
                                zorder=-1)
                    ax.axvline(real_values.mean(), ls='--', color="#808080",
                               zorder=-1)
                except RuntimeError:
                    pass

            if max(predict_values) > 0:
                try:
                    sns.kdeplot(predict_values, shade=True, color=color, ax=ax,
                                zorder=1)
                    ax.axvline(predict_values.mean(), ls='--', color=color,
                               zorder=1)
                except RuntimeError:
                    pass

            ax.set_title(title, fontsize=20, fontweight='bold')
            ax.set_xlabel("weights", fontsize=16, fontweight='bold')
            ax.tick_params(labelsize=14)

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "weights_distributions.{}".format(name, self.extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory path: {}".format(self.input_path))
        print("  > Input file suffix: {}".format(self.input_suffix))
        print("  > Reference profile path: {}".format(self.ref_profile_path))
        print("  > Gene info path: {}".format(self.gene_info_path))
        print("  > Combine: {}".format(self.combine))
        print("  > Output directory path: {}".format(self.outdir))
        print("")
