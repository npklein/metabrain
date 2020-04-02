"""
File:         workers2.py
Created:      2020/04/01
Last Changed: 2020/04/02
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
                   start_q, job_q, result_q, sleep_time, max_end_time, verbose):
    start_time = time.time()
    start_time_str = datetime.fromtimestamp(start_time).strftime(
        "%d-%m-%Y, %H:%M:%S")
    print("[worker {:2d}]\tstarted [{}].".format(worker_id, start_time_str))

    # Wait until sample orders are received.
    sample_orders = None
    while sample_orders is None:
        if max_end_time <= time.time():
            break

        try:
            if start_q.empty():
                break
            sample_orders = start_q.get(True, timeout=1)
        except queue.Empty:
            break

    if sample_orders is None:
        return

    if verbose:
        print("[worker {:2d}]\treceived sample orders.".format(worker_id))

    # Load the covariate file.
    tech_cov_df, cov_df = load_covariate_file(cov_inpath, tech_covs)

    if verbose:
        print("[worker {:2d}]\tloaded covariate file.".format(worker_id))
        print("[worker {:2d}]\t\t"
              "technical covariates: {}".format(worker_id, tech_cov_df.shape))
        print("[worker {:2d}]\t\t"
              "relevant covariates:  {}".format(worker_id, cov_df.shape))

    # Push the header to the manager. This is also the message that confirms
    # the worker is ready to start.
    result_q.put((worker_id, "online", -1, None, ["-"] + list(cov_df.index)))

    # Start infinite loop.
    while True:
        if max_end_time <= time.time():
            break

        try:
            if job_q.empty():
                time.sleep(sleep_time)

            # Get my work.
            schedule = job_q.get(True, timeout=1)
            if worker_id not in schedule.keys() or schedule[worker_id] is None:
                result_q.put((worker_id, "waiting", None, None, None))
                time.sleep(sleep_time)
                continue
            my_work = schedule[worker_id]
            if verbose:
                print("[worker {:2d}]\tFound work".format(worker_id))

            # Combine start indices.
            start_index = min(my_work.keys())
            nrows = max(my_work.keys()) - start_index + 1

            if verbose:
                print("[worker {:2d}]\tloading data {}-{}".format(worker_id,
                                                                 start_index,
                                                                 start_index + nrows))

            # Load the genotype and expression data.
            try:
                geno_df = pd.read_csv(geno_inpath,
                                      sep="\t",
                                      header=0,
                                      index_col=0,
                                      skiprows=[i for i in
                                                range(1, start_index)],
                                      nrows=nrows)

                expr_df = pd.read_csv(expr_inpath,
                                      sep="\t",
                                      header=0,
                                      index_col=0,
                                      skiprows=[i for i in
                                                range(1, start_index)],
                                      nrows=nrows)
            except pandas.io.common.EmptyDataError:
                break

            # Check if the covariate data still exists. If not, reload it.
            try:
                _ = tech_cov_df.shape
                _ = cov_df.shape
            except:
                print("[worker {:2d}]\tcovariate table was lost. "
                      "Reloading data..".format(worker_id))
                tech_cov_df, cov_df = load_covariate_file(cov_inpath, tech_covs)

            # Replace -1 with NaN in the genotype dataframe.
            geno_df.replace(-1, np.nan, inplace=True)

            # loop over eQTLs.
            for i, eqtl_index in enumerate(my_work.keys()):
                if verbose:
                    print("[worker {:2d}]\tworking on 'eQTL_{}'.".format(worker_id, eqtl_index))

                # Get the missing genotype indices.
                indices = np.arange(geno_df.shape[1])
                eqtl_indices = indices[~geno_df.iloc[i, :].isnull().values]
                n = len(eqtl_indices)

                # Subset the row and present samples for this eQTL.
                genotype = geno_df.iloc[i, eqtl_indices].copy()
                genotype_all = geno_df.iloc[i, :].copy()
                expression = expr_df.iloc[i, eqtl_indices].copy()
                technical_covs = tech_cov_df.iloc[:, eqtl_indices].copy()

                # Check if SNP index are identical.
                if genotype.name != expression.name:
                    print("[worker {:2d}]\tindices do not match.".format(worker_id))
                    continue

                # Create the null model.
                null_matrix = technical_covs.mul(genotype, axis=1)
                df_null, rss_null = create_model(null_matrix.T, expression)

                # Loop over each sample order.
                for order_id, sample_order in enumerate(sample_orders):
                    if order_id not in my_work[eqtl_index]:
                        continue
                    if verbose:
                        print("[worker {:2d}]\tworking on 'sample order {}'.".format(worker_id, order_id))

                    # Loop over the covariates.
                    snp_pvalues = []
                    for j in range(cov_df.shape[0]):
                        cov_name = cov_df.index[j]

                        # if verbose:
                        #     print("[worker {:2d}]\tworking on covariate "
                        #           "{}/{}".format(worker_id, j + 1,
                        #                          cov_df.shape[0]))

                        # Calculate the interaction of the covariate of
                        # interest.
                        covariate_all = cov_df.iloc[j, :].copy()
                        interaction_all = covariate_all * genotype_all
                        interaction_all.name = cov_name

                        # Reorder the series.
                        interaction_all = interaction_all.reindex(interaction_all.index[sample_order])

                        # Remove missing values.
                        interaction_all = interaction_all.dropna()

                        # Set the identical names as the null matrix.
                        interaction_all.index = null_matrix.columns

                        # Create the alternative matrix.
                        alt_matrix = null_matrix.copy()
                        alt_matrix = alt_matrix.append(interaction_all)

                        # Create the alternative model.
                        df_alt, rss_alt = create_model(alt_matrix.T, expression)

                        # Compare the models.
                        pvalue = compare_models(rss_null, rss_alt,
                                                df_null, df_alt,
                                                n)

                        # Safe.
                        snp_pvalues.append(pvalue)

                    # Push the results.
                    if verbose:
                        print("[worker {:2d}]\tpushing results".format(worker_id))
                    result_q.put((worker_id, "result", eqtl_index, order_id, [genotype.name] + snp_pvalues))

            # Wait before getting a new job.
            time.sleep(sleep_time * 2)

        except queue.Empty:
            break

    # Calculate the worker process time.
    end_time = time.time()
    end_time_str = datetime.fromtimestamp(end_time).strftime(
        "%d-%m-%Y, %H:%M:%S")
    run_time_min, run_time_sec = divmod(end_time - start_time, 60)
    run_time_hour, run_time_min = divmod(run_time_min, 60)
    print("[worker {:2d}]\tstopped [{}], worked for {} hour(s), "
          "{} minute(s), and {} second(s).".format(worker_id, end_time_str,
                                                   int(run_time_hour),
                                                   int(run_time_min),
                                                   int(run_time_sec)))


def load_covariate_file(cov_inpath, tech_covs):
    cov_df = pd.read_csv(cov_inpath,
                         sep="\t",
                         header=0,
                         index_col=0)

    tech_cov_df = cov_df.loc[tech_covs, :]
    cov_df.drop(tech_covs, inplace=True)

    return tech_cov_df, cov_df


def create_model(X, y):
    regressor = LinearRegression(n_jobs=1)
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
    # stats.f.sf(1827.95, dfn=1, dfd=3661) = 5e-324
    # stats.f.sf(9.9e-12, dfn=1, dfd=3661) = 0.9999974567714613

    # stats.f.cdf(69, dfn=1, dfd=3661) = 0.9999999999999999
    # stats.f.cdf(1e-320, dfn=1, dfd=3661) = 1.0730071046473278e-160

    return stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2))
