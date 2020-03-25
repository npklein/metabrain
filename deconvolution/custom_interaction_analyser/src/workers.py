"""
File:         workers.py
Created:      2020/03/24
Last Changed: 2020/03/25
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
    print("Worker: {} started".format(worker_id))

    # Load the covariate file.
    cov_df = pd.read_csv(cov_inpath,
                         sep="\t",
                         header=0,
                         index_col=0)

    tech_cov = cov_df.loc[tech_covs, :]
    cov_df.drop(tech_covs, inplace=True)

    output_q.put((None, [-1] + ["-"] + list(cov_df.index)))

    while True:
        if max_end_time <= time.time():
            break

        try:
            if input_q.empty():
                break
            start_index = input_q.get(True, timeout=1)

            try:
                geno_df = pd.read_csv(geno_inpath,
                                      sep="\t",
                                      header=0,
                                      index_col=0,
                                      skiprows=[i for i in range(1, start_index)],
                                      nrows=chunk_size)

                expr_df = pd.read_csv(expr_inpath,
                                      sep="\t",
                                      header=0,
                                      index_col=0,
                                      skiprows=[i for i in range(1, start_index)],
                                      nrows=chunk_size)

            except pandas.io.common.EmptyDataError:
                break

            # Replace -1 with NaN.
            geno_df.replace(-1, np.nan, inplace=True)

            # loop over eQTLs.
            for i in range(chunk_size):
                print("Worker: {}, working on eQTL {}/{}".format(worker_id,
                                                                 start_index + i,
                                                                 (start_index + chunk_size) - 1))
                # Get the missing genotype indices.
                indices = np.arange(geno_df.shape[1])
                sample_mask = indices[~geno_df.iloc[i, :].isnull().values]
                n = len(sample_mask)

                # Subset the row and present samples for this eQTL.
                genotype = geno_df.iloc[i, sample_mask].copy()
                expression = expr_df.iloc[i, sample_mask].copy()
                technical_covs = tech_cov.iloc[:, sample_mask].copy()
                covariates = cov_df.iloc[:, sample_mask].copy()

                # Check if SNP index are identical.
                if genotype.name != expression.name:
                    print("Indices do not match.")
                    continue

                # Create the null model.
                null_matrix = technical_covs.mul(genotype, axis=1)
                df_null, rss_null = create_model(null_matrix.T, expression)

                # Loop over the covariates.
                snp_zscores = []
                snp_adj_zscores = []
                for j in range(covariates.shape[0]):
                    print("Worker: {}, working on covariate {}/{}".format(worker_id,
                                                                          j + 1,
                                                                          covariates.shape[0]))
                    # Append to covariate.
                    name = covariates.index[j]
                    cov_of_interest = covariates.iloc[j, :]
                    alt_matrix = null_matrix.copy()
                    alt_matrix.loc[name, :] = cov_of_interest * genotype

                    # Create the alternative model.
                    df_alt, rss_alt = create_model(alt_matrix.T, expression)

                    # Compare the models.
                    pvalue, zscore = compare_models(rss_null, rss_alt, df_null, df_alt, n)

                    # Safe.
                    snp_zscores.append(zscore)

                    # Permutate.
                    if n_permutations > 0:
                        perm_count = 0
                        for k in range(n_permutations):
                            # Shuffle the genotypes.
                            inter_array = cov_of_interest * genotype
                            index = inter_array.index
                            inter_array = inter_array.sample(frac=1)
                            inter_array.index = index

                            # Append to covariate.
                            perm_matrix = null_matrix.copy()
                            perm_matrix.loc[name, :] = inter_array

                            # Create the alternative model.
                            df_alt, rss_alt = create_model(perm_matrix.T, expression)

                            # Compare the models.
                            perm_pvalue, _ = compare_models(rss_null, rss_alt, df_null, df_alt, n)
                            if pvalue < perm_pvalue:
                                perm_count += 1

                        # Calculate the emperical FDR.
                        adj_pvalue = perm_count / n_permutations
                        adj_zscore = get_z_score(adj_pvalue)

                        # Safe.
                        snp_adj_zscores.append(adj_zscore)

                # Push the results.
                print("Worker: {} pushing results".format(worker_id))
                output_q.put(("default", [start_index + i] + [genotype.name] + snp_zscores))
                output_q.put(("adjusted", [start_index + i] + [genotype.name] + snp_adj_zscores))

        except queue.Empty:
            break

    print("Worker: {} stopped".format(worker_id))


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
    z_value = get_z_score(p_value)

    return p_value, z_value


def calc_f_value(rss1, rss2, df1, df2, n):
    return((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))


def get_p_value(f_value, df1, df2, n):
    p_value = 1 - stats.f.cdf(f_value,
                              dfn=(df2 - df1),
                              dfd=(n - df2))
    if p_value < 2.2250738585072014e-308:
        p_value = 2.2250738585072014e-308

    return p_value


def get_z_score(p_value):
    return stats.norm.ppf((1 - p_value))