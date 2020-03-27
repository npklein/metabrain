"""
File:         workers.py
Created:      2020/03/24
Last Changed: 2020/03/27
Author(s):    M.Vochteloo

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
from datetime import datetime
import queue
import time

# Third party imports.
import numpy as np
import pandas as pd
import pandas.io.common
from sklearn.linear_model import LinearRegression
import scipy.stats as stats


# Local application imports.


def process_worker(worker_id, cov_inpath, geno_inpath, expr_inpath, tech_covs,
                   input_q, output_q, max_end_time, chunk_size, n_permutations):
    start_time = time.time()
    start_time_str = datetime.fromtimestamp(start_time).strftime(
        "%d-%m-%Y, %H:%M:%S")
    print("Worker: {} started [{}]".format(worker_id, start_time_str))

    # Create time lists.
    eqtl_process_times = []

    # Load the covariate file.
    cov_df = pd.read_csv(cov_inpath,
                         sep="\t",
                         header=0,
                         index_col=0)

    tech_cov = cov_df.loc[tech_covs, :]
    cov_df.drop(tech_covs, inplace=True)
    cov_columns = cov_df.columns

    output_q.put((None, [-1] + ["-"] + list(cov_df.index)))

    while True:
        if max_end_time <= time.time():
            break

        try:
            if input_q.empty():
                break
            (start_index, shuffle_type, column_indices) = input_q.get(True, timeout=1)
            # print("Worker: {} start_index: {}\tshuffle type: {}\t"
            #       "column_indices: {}".format(worker_id, start_index,
            #                                   shuffle_type, column_indices[:5]))

            try:
                geno_df = pd.read_csv(geno_inpath,
                                      sep="\t",
                                      header=0,
                                      index_col=0,
                                      skiprows=[i for i in
                                                range(1, start_index)],
                                      nrows=chunk_size)

                expr_df = pd.read_csv(expr_inpath,
                                      sep="\t",
                                      header=0,
                                      index_col=0,
                                      skiprows=[i for i in
                                                range(1, start_index)],
                                      nrows=chunk_size)

            except pandas.io.common.EmptyDataError:
                break

            # Reorder the genotype dataframe.
            geno_columns = geno_df.columns
            geno_df = geno_df.iloc[:, column_indices]
            geno_df.columns = geno_columns

            # Reorder the covariate dataframe.
            tmp_cov_df = cov_df.copy()
            tmp_cov_df = tmp_cov_df.iloc[:, column_indices]
            tmp_cov_df.columns = cov_columns

            # Replace -1 with NaN in the genotype dataframe.
            geno_df.replace(-1, np.nan, inplace=True)

            # loop over eQTLs.
            for i in range(chunk_size):
                eqtl_start_time = time.time()
                # print("Worker: {}, working on eQTL {}/{}".format(worker_id,
                #                                                  start_index + i,
                #                                                  (start_index + chunk_size) - 1))
                # Get the missing genotype indices.
                indices = np.arange(geno_df.shape[1])
                sample_mask = indices[~geno_df.iloc[i, :].isnull().values]
                n = len(sample_mask)

                # Subset the row and present samples for this eQTL.
                genotype = geno_df.iloc[i, sample_mask].copy()
                expression = expr_df.iloc[i, sample_mask].copy()
                technical_covs = tech_cov.iloc[:, sample_mask].copy()
                covariates = tmp_cov_df.iloc[:, sample_mask].copy()

                # Check if SNP index are identical.
                if genotype.name != expression.name:
                    print("Indices do not match.")
                    continue

                # Create the null model.
                null_matrix = technical_covs.mul(genotype, axis=1)
                df_null, rss_null = create_model(null_matrix.T, expression)

                # Loop over the covariates.
                snp_pvalues = []
                for j in range(covariates.shape[0]):
                    # print("Worker: {}, working on covariate {}/{}".format(
                    #     worker_id,
                    #     j + 1,
                    #     covariates.shape[0]))
                    # Append to covariate.
                    cov_name = covariates.index[j]
                    cov_of_interest = covariates.iloc[j, :]
                    alt_matrix = null_matrix.copy()
                    alt_matrix.loc[cov_name, :] = cov_of_interest * genotype

                    # Create the alternative model.
                    df_alt, rss_alt = create_model(alt_matrix.T, expression)

                    # Compare the models.
                    pvalue = compare_models(rss_null, rss_alt, df_null, df_alt,
                                            n)

                    # Safe.
                    snp_pvalues.append(pvalue)

                # Push the results.
                # print("Worker: {} pushing results".format(worker_id))
                output_q.put((shuffle_type, [start_index + i] + [genotype.name] + snp_pvalues))

                # Safe the time.
                eqtl_process_time = time.time() - eqtl_start_time
                eqtl_process_times.append(eqtl_process_time)

        except queue.Empty:
            break

    # Print the process times.
    print("eqtl_process_times <- c({})".format(
        ', '.join([str(round(x, 2)) for x in eqtl_process_times])))

    # Calculate the worker process time.
    end_time = time.time()
    end_time_str = datetime.fromtimestamp(end_time).strftime(
        "%d-%m-%Y, %H:%M:%S")
    run_time_min, run_time_sec = divmod(end_time - start_time, 60)
    run_time_hour, run_time_min = divmod(run_time_min, 60)
    print("Worker: {} stopped [{}], ran for {} hours, "
          "{} minutes and {} seconds.".format(worker_id, end_time_str,
                                              int(run_time_hour),
                                              int(run_time_min),
                                              int(run_time_sec)))


def perform_permutations(null_matrix, cov_of_interest, genotype, expression,
                         cov_name, rss_null, df_null, pvalue, n, n_permutations):
    p_values = []
    pvalue_rank = 0
    for k in range(n_permutations):
        # Shuffle the genotypes.
        inter_array = cov_of_interest * genotype
        index = inter_array.index
        inter_array = inter_array.sample(frac=1)
        inter_array.index = index

        # Append to covariate.
        perm_matrix = null_matrix.copy()
        perm_matrix.loc[cov_name, :] = inter_array

        # Create the alternative model.
        df_alt, rss_alt = create_model(perm_matrix.T, expression)

        # Compare the models.
        perm_pvalue, _ = compare_models(rss_null, rss_alt, df_null, df_alt,
                                        n)
        if perm_pvalue < pvalue:
            pvalue_rank += 1
        p_values.append(perm_pvalue)

    # Calculate the emperical FDR.
    adj_pvalue = calc_adj_p_value(pvalue_rank, n_permutations)

    return adj_pvalue


def create_model(X, y):
    regressor = LinearRegression()
    regressor.fit(X, y)
    y_hat = regressor.predict(X)

    # Calculate the statistics of the model.
    degrees_freedom = len(X.columns) + 1
    residual_squared_sum = ((y - y_hat) ** 2).sum()

    return degrees_freedom, residual_squared_sum


def compare_models(rss1, rss2, df1, df2, n):
    f_value = calc_f_value(rss1, rss2, df1, df2, n)
    p_value = get_p_value(f_value, df1, df2, n)

    return p_value


def calc_f_value(rss1, rss2, df1, df2, n):
    if df1 >= df2:
        return np.nan
    if df2 >= n:
        return np.nan

    return ((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))


def get_p_value(f_value, df1, df2, n):
    p_value = 1 - stats.f.cdf(f_value,
                              dfn=(df2 - df1),
                              dfd=(n - df2))
    if p_value < 1e-16:
        p_value = 1e-16

    return p_value


def calc_adj_p_value(rank, n_permutations):
    adj_pvalue = (rank + 1) / (n_permutations + 1)
    if adj_pvalue >= 1:
        adj_pvalue = 1 - 1e-16
    elif adj_pvalue < 1e-16:
        adj_pvalue = 1e-16
    return adj_pvalue
