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
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy.optimize import nnls
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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

        self.palette = {
            "NEU": "#0072B2",
            "OLI": "#009E73",
            "END": "#CC79A7",
            "MIC": "#E69F00",
            "AST": "#D55E00"
        }

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        ref_profile_df = self.load_reference_profile()
        ref_profile_df = self.filter(ref_profile_df, 5)
        ref_profile_df = self.perform_log2_transform(ref_profile_df)
        if ref_profile_df.values.min() < 0:
            ref_profile_df = self.perform_shift(ref_profile_df)

        cell_types = list(ref_profile_df.columns)
        input_cell_types = list(ref_profile_df.columns)
        if self.combine:
            combine_list = [x.split("-") for x in self.combine]
            for old1, old2, new in combine_list:
                input_cell_types.remove(new)
                input_cell_types.append(old1)
                input_cell_types.append(old2)

        print("")
        print("### Step2 ###")
        samples, matrices = self.load_ct_matrices(cell_types=input_cell_types,
                                                  genes=ref_profile_df.index)

        print("")
        print("### Step3 ###")
        cc_df = self.load_cell_counts(cell_types=input_cell_types)
        print(cc_df)

        print("")
        print("### Step4 ###")
        sample_matrices = self.merge_samples(samples=samples,
                                             cell_types=cell_types,
                                             matrices=matrices)

        print("")
        print("### Step5 ###")
        real_data = []
        predict_data = []
        diff_data = []
        for sample in samples:
            if sample not in cc_df:
                continue
            real_weights = cc_df.loc[:, sample]

            _, bulk_expression = self.create_artifical_bulk_data(df=sample_matrices[sample],
                                                                 weights=real_weights)

            decon_weights = self.deconvolute(ref_profile_df, bulk_expression)

            # weights_df = pd.DataFrame({"real": real_weights,
            #                            "predict": decon_weights})
            # weights_df.index = ref_profile_df.columns
            # weights_df["diff"] = (weights_df["real"] - weights_df["predict"]).abs()
            # print(weights_df)

            real_data.append(real_weights)
            predict_data.append(decon_weights)
            diff_data.append([a-b for a, b in zip(real_weights, decon_weights)])

        real_weights_df = pd.DataFrame(real_data, columns=cell_types)
        decon_weights_df = pd.DataFrame(predict_data, columns=cell_types)
        diff_df = pd.DataFrame(diff_data, columns=cell_types).abs()

        print("")
        print("### Step6 ###")
        print("Visualizing differences")
        real_weights_df_m = real_weights_df.melt()
        real_weights_df_m["hue"] = "real"
        decon_weights_df_m = decon_weights_df.melt()
        decon_weights_df_m["hue"] = "predict"

        dfm = pd.concat([real_weights_df_m, decon_weights_df_m], axis=0)
        self.plot_distributions(dfm, name="weights_comparison")

        print("")
        print("### Results ###")
        print("Difference between real and predicted weights:")
        print(diff_df.describe())

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

    def load_ct_matrices(self, cell_types, genes):
        print("Loading cell type matrices")
        columns = None
        matrices = {}

        for cell_type in cell_types:
            inpath = os.path.join(self.input_path, cell_type + self.input_suffix)
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

        if self.combine is not None:
            for pair in self.combine:
                left_ct, right_ct, name = pair.split("-")

                left_df = full_matrices[left_ct]
                right_df = full_matrices[right_ct]

                full_matrices[name] = left_df.add(right_df, fill_value=0).div(2)

        return columns, full_matrices

    def load_cell_counts(self, cell_types):
        df = load_dataframe(self.cell_counts_path,
                            sep="\t",
                            header=0,
                            index_col=0)

        df = df.loc[cell_types, :]
        print(df)

        if self.combine is not None:
            for pair in self.combine:
                left_ct, right_ct, name = pair.split("-")

                df.loc[name, :] = df.loc[left_ct, :] + df.loc[right_ct, :]
                df.drop(left_ct, axis=0, inplace=True)
                df.drop(right_ct, axis=0, inplace=True)

        return df / df.sum(axis=0)

    @staticmethod
    def merge_samples(samples, cell_types, matrices):
        sample_matrices = {}

        for sample in samples:
            print("\tCombining cell type expression for sample {}".format(sample))
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
        return weights / np.sum(weights)

    def plot_distributions(self, df, name=""):
        sns.set_style("ticks")
        fig, axes = plt.subplots(nrows=1,
                                 ncols=len(df["variable"].unique()),
                                 figsize=(8*len(df["variable"].unique()), 6))

        for i, ct in enumerate(df["variable"].unique()):
            ax = axes[i]
            sns.despine(fig=fig, ax=ax)

            real_values = df.loc[(df["variable"] == ct) & (df["hue"] == "real"), "value"]
            real_values.name = "real"
            predict_values = df.loc[(df["variable"] == ct) & (df["hue"] == "predict"), "value"]
            predict_values.name = "predict"

            sns.kdeplot(real_values, shade=True, color="#808080", ax=ax,
                        zorder=-1)
            ax.axvline(real_values.mean(), ls='--', color="#808080",
                       zorder=-1)

            sns.kdeplot(predict_values, shade=True, color=self.palette[ct], ax=ax,
                        zorder=1)
            ax.axvline(predict_values.mean(), ls='--', color=self.palette[ct],
                       zorder=1)

            ax.set_title(ct, fontsize=20, fontweight='bold')
            ax.set_xlabel("weights", fontsize=16, fontweight='bold')
            ax.tick_params(labelsize=14)

        plt.tight_layout()
        fig.savefig(os.path.join(self.outdir, "{}_distributions.{}".format(name, self.extension)))
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
