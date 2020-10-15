"""
File:         combiner.py
Created:      2020/10/15
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
from itertools import groupby, count
import pickle
import glob
import time
import os

# Third party imports.
import numpy as np
import pandas as pd
from scipy import stats

# Local application imports.
from local_settings import LocalSettings
from utilities import save_dataframe


class Combine:
    def __init__(self, name, settings_file):
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        settings = LocalSettings(current_dir, settings_file)

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir, name)

        # Get the needed settings.
        self.cov_outdir = settings.get_setting("covariates_folder")
        self.tech_cov_outdir = settings.get_setting("technical_covariates_folder")
        self.pvalues_filename = settings.get_setting("real_pvalues_pickle_filename")
        self.snp_coef_filename = settings.get_setting("snp_coef_pickle_filename")
        self.snp_std_err_filename = settings.get_setting("snp_std_err_pickle_filename")
        self.snp_tvalue_filename = settings.get_setting("snp_tvalue_pickle_filename")
        self.inter_coef_filename = settings.get_setting("inter_coef_pickle_filename")
        self.inter_std_err_filename = settings.get_setting("inter_std_err_pickle_filename")
        self.inter_tvalue_filename = settings.get_setting("inter_tvalue_pickle_filename")
        self.perm_pvalues_filename = settings.get_setting("permuted_pvalues_pickle_filename")

    def start(self):
        print("Starting interaction analyser - combine and plot.")
        self.print_arguments()

        # Start the timer.
        start_time = time.time()

        for outdir in [self.cov_outdir, self.tech_cov_outdir]:
            full_outdir = os.path.join(self.outdir, outdir)
            if not os.path.exists(full_outdir):
                print("Error, output directory does not exist.")

            self.work(full_outdir)

        # Print the time.
        run_time_min, run_time_sec = divmod(time.time() - start_time, 60)
        run_time_hour, run_time_min = divmod(run_time_min, 60)
        print("finished in  {} hour(s), {} minute(s) and "
              "{} second(s).".format(int(run_time_hour),
                                     int(run_time_min),
                                     int(run_time_sec)), flush=True)

    def work(self, workdir):
        dataframes = {}

        print("### Step 1 ###", flush=True)
        print("Combine pickle files into dataframe.")
        for filename in [self.pvalues_filename,
                         self.snp_coef_filename,
                         self.snp_std_err_filename,
                         self.inter_coef_filename,
                         self.inter_std_err_filename]:

            print("\t Loading {} data.".format(filename), flush=True)
            columns, data = self.combine_pickles(workdir,
                                                 filename,
                                                 columns=True)
            print(len(data))

            print("Creating {} dataframe.".format(filename), flush=True)
            df = self.create_df(data, columns)

            print("Saving {} dataframe.".format(filename), flush=True)
            save_dataframe(df=df,
                           outpath=os.path.join(workdir,
                                                "{}_table.txt.gz".format(filename)),
                           header=True, index=True)

            dataframes[filename] = df.copy()

            del columns, data, df

        print("")
        print("### Step 2 ###", flush=True)
        print("Calculate t-values")
        for coef_file, std_err_file, tvalue_file in [[self.snp_coef_filename, self.snp_std_err_filename, self.snp_tvalue_filename],
                                                     [self.inter_coef_filename, self.inter_std_err_filename, self.inter_tvalue_filename]]:
            coef_df = dataframes[coef_file]
            std_err_df = dataframes[std_err_file]
            tvalue_df = coef_df / std_err_df

            print("Saving {} dataframe.".format(tvalue_file), flush=True)
            save_dataframe(df=tvalue_df,
                           outpath=os.path.join(workdir,
                                                "{}_table.txt.gz".format(tvalue_file)),
                           header=True, index=True)

        print("")
        print("### Step 3 ###", flush=True)
        print("Calculate z-scores")

        pvalue_df = dataframes[self.pvalues_filename]

        print("Creating z-score dataframe.", flush=True)
        zscore_df = self.create_zscore_df(pvalue_df.copy())

        print("Saving z-score dataframe.", flush=True)
        save_dataframe(df=zscore_df,
                       outpath=os.path.join(workdir,
                                            "zscore_table.txt.gz"),
                       header=True, index=True)

        print("")
        print("### Step 4 ###", flush=True)
        print("Calculate FDR-values")

        print("Loading permutation pvalue data.", flush=True)
        _, perm_pvalues = self.combine_pickles(workdir,
                                               self.perm_pvalues_filename)
        n_permutations = 0
        fdr_df = None
        if len(perm_pvalues) > 0:
            n_permutations = len(perm_pvalues) / len(pvalue_df.index)

            # TODO

        print("")
        print("Program completed.", flush=True)

    @staticmethod
    def combine_pickles(indir, filename, columns=False):
        # Declare variables.
        col_list = None
        data = []

        # Combine the found files.
        for i, fpath in enumerate(glob.glob(os.path.join(indir, filename,
                                                         filename + "*.pkl"))):
            with open(fpath, "rb") as f:
                try:
                    content = pickle.load(f)
                    if columns and i == 0:
                        col_list = content[0]
                    data.extend(content[1:])
                except EOFError:
                    print("\tEOFError in: {} ".format(os.path.basename(fpath)))
            f.close()

        return col_list, data

    def create_df(self, data, columns):
        # Crate the data frame.
        df = pd.DataFrame(data, columns=columns)
        order_col = df.columns[0]
        df.sort_values(by=order_col, inplace=True)
        print("\tInput shape: {}".format(df.shape))

        # Check what eQTLs are present, missing, duplicated.
        reference = set(np.arange(df[order_col].min(), df[order_col].max() + 1))
        present = set(df[order_col])
        print("\tPresent indices: {}".format(self.group_consecutive_numbers(present)))
        missing = list(set(df[order_col]).symmetric_difference(reference))
        print("\tMissing indices: {}".format(self.group_consecutive_numbers(missing)))
        duplicated = list(df.loc[df.duplicated([order_col]), order_col].values)
        print("\tDuplicate indices: {}".format(self.group_consecutive_numbers(duplicated)))

        # Remove duplicate eQTLs.
        df.drop_duplicates(subset=order_col, keep="first", inplace=True)

        # Insert missing eQTLs.
        if len(missing) > 0:
            print("\tInserting missing indices")
            missing_df = pd.Series(missing).to_frame()
            missing_df.columns = [order_col]
            df = pd.concat([df, missing_df], axis=0).reset_index(drop=True)
            df.sort_values(by=order_col, inplace=True)

        # Remove the eQTL order column and set the SNPName as index.
        df.drop([order_col], axis=1, inplace=True)
        df.set_index("-", inplace=True)
        df = df.T
        print("\tOutput shape: {}".format(df.shape))

        return df

    @staticmethod
    def group_consecutive_numbers(numbers):
        groups = groupby(numbers, key=lambda item, c=count(): item - next(c))
        tmp = [list(g) for k, g in groups]
        return [str(x[0]) if len(x) == 1 else "{}-{}".format(x[0], x[-1]) for x in tmp]

    def create_zscore_df(self, df):
        count = 1
        total = df.shape[0] * df.shape[1]

        zscore_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                zscore = self.get_z_score(pvalue)
                zscore_df.iloc[row_index, col_index] = zscore
                count += 1

        return zscore_df

    @staticmethod
    def get_z_score(p_value):
        if p_value > (1.0 - 1e-16):
            p_value = (1.0 - 1e-16)
        if p_value < 1e-323:
            p_value = 1e-323

        # The lower and upper limit of stats.norm.sf
        # stats.norm.isf((1 - 1e-16)) = -8.209536151601387
        # stats.norm.isf(1e-323) = 38.44939448087599
        return stats.norm.isf(p_value)

    @staticmethod
    def create_perm_fdr_df(df, pvalues, perm_pvalues, n_perm):
        count = 1
        total = df.shape[0] * df.shape[1]

        max_signif_pvalue = -np.inf

        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue)
                perm_rank = bisect_left(perm_pvalues, pvalue)
                if (rank > 0) and (perm_rank > 0):
                    fdr_value = (perm_rank / n_perm) / rank
                    if fdr_value > 1:
                        fdr_value = 1
                else:
                    fdr_value = 0
                fdr_df.iloc[row_index, col_index] = fdr_value

                if fdr_value < 0.05 and pvalue > max_signif_pvalue:
                    max_signif_pvalue = pvalue

                count += 1

        return fdr_df, max_signif_pvalue

    @staticmethod
    def create_bh_fdr_df(df, pvalues):
        count = 1
        total = df.shape[0] * df.shape[1]

        max_signif_pvalue = -np.inf

        m = np.count_nonzero(~np.isnan(pvalues))
        fdr_df = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        fdr_values = []
        for row_index in range(df.shape[0]):
            for col_index in range(df.shape[1]):
                if (count == 1) or (count % int(total / 10) == 0) or (count == total):
                    print("\tworking on {}/{} [{:.2f}%]".format(count,
                                                                total,
                                                                (100 / total) * count))
                pvalue = df.iloc[row_index, col_index]
                rank = bisect_left(pvalues, pvalue) + 1
                fdr_value = pvalue * (m / rank)
                if fdr_value > 1:
                    fdr_value = 1
                fdr_values.append((rank, row_index, col_index, fdr_value))
                fdr_df.iloc[row_index, col_index] = fdr_value

                if fdr_value < 0.05 and pvalue > max_signif_pvalue:
                    max_signif_pvalue = pvalue

                count += 1

        # Make sure the BH FDR is a monotome function. This goes through
        # the FDR values backwords and make sure that the next FDR is always
        # lower or equal to the previous.
        fdr_values = sorted(fdr_values, key=lambda x: x[0], reverse=True)
        prev_fdr_value = None
        for (rank, row_index, col_index, fdr_value) in fdr_values:
            if prev_fdr_value is not None and fdr_value > prev_fdr_value:
                fdr_df.iloc[row_index, col_index] = prev_fdr_value
                prev_fdr_value = prev_fdr_value
            else:
                prev_fdr_value = fdr_value

        return fdr_df, max_signif_pvalue



    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("")
