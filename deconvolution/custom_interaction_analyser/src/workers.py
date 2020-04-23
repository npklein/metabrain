"""
File:         workers.py
Created:      2020/04/01
Last Changed: 2020/04/23
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
import scipy.stats as stats
import statsmodels.api as sm
from functools import reduce

# Local application imports.


def process_worker(worker_id, cov_inpath, geno_inpath, expr_inpath, tech_covs,
                   start_q, job_q, result_q, sleep_time, max_time_unresponsive,
                   max_end_time, verbose):
    start_time_str = datetime.fromtimestamp(time.time()).strftime(
        "%d-%m-%Y, %H:%M:%S")
    print("[worker {:2d}]\tstarted [{}]".format(worker_id, start_time_str),
          flush=True)

    # Define the update frequency.
    update_frequency = int(max_time_unresponsive / 4.1)

    # First receive the permutation order from the start queue.
    sample_orders = None
    while sample_orders is None:
        # Break the loop of the maximum runtime has been reached. This prevents
        # the workers to continue indefinitely.
        if max_end_time <= time.time():
            break

        try:
            if start_q.empty():
                break
            sample_orders = start_q.get(True, timeout=1)
        except queue.Empty:
            break

    # The program cannot start if it hasn't gotten the sample orders.
    if sample_orders is None:
        return

    if verbose:
        print("[worker {:2d}]\treceived sample orders".format(worker_id),
              flush=True)

    # Load the covariate file. Split the matrix into technical covariates and
    # covariates of interest.
    tech_cov_df, cov_df = load_covariate_file(cov_inpath, tech_covs)

    if verbose:
        print("[worker {:2d}]\tloaded covariate file".format(worker_id),
              flush=True)
        print("[worker {:2d}]\t\t"
              "technical covariates: {}".format(worker_id, tech_cov_df.shape),
              flush=True)
        print("[worker {:2d}]\t\t"
              "relevant covariates:  {}".format(worker_id, cov_df.shape),
              flush=True)

    # Push the header to the manager. This is also the message that confirms
    # the worker is ready to start.
    result_q.put((worker_id, "online", -1, None, ["-"] + list(cov_df.index)))
    last_message = int(time.time())

    # Start infinite loop.
    while True:
        # Break the loop of the maximum runtime has been reached. This prevents
        # the workers to continue indefinitely.
        if max_end_time <= time.time():
            break

        try:
            if job_q.empty():
                time.sleep(sleep_time)

            # Get the schedule from the job queue.
            schedule = job_q.get(True, timeout=1)

            # Check if there is work for me, if not, send a waiting result
            # back to the manager
            if worker_id not in schedule.keys() or schedule[worker_id] is None:
                now = int(time.time())
                if (now - last_message) > update_frequency:
                    if verbose:
                        print("[worker {:2d}]\tSending "
                              "heartbeat".format(worker_id), flush=True)
                    result_q.put((worker_id, "waiting", None, None, None))
                    last_message = now
                time.sleep(sleep_time)
                continue
            # The work has format: {eqtl_id [sample order id, ...], ...}
            my_work = schedule[worker_id]
            if verbose:
                print("[worker {:2d}]\tFound work".format(worker_id),
                      flush=True)

            # Find the start index for reading the input files. The lowest
            # key is the start index for reading the input file. The end nrows
            # is the range of start indices + 1.
            start_index = min(my_work.keys())
            nrows = max(my_work.keys()) - start_index + 1

            if verbose:
                print("[worker {:2d}]\tloading data "
                      "{}-{}".format(worker_id, start_index,
                                     start_index + nrows - 1), flush=True)

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

            # Replace -1 with NaN in the genotype dataframe. This way we can
            # drop missing values.
            geno_df.replace(-1, np.nan, inplace=True)

            # loop over the eQTLs.
            for row_index, eqtl_index in enumerate(my_work.keys()):
                if verbose:
                    print("[worker {:2d}]\tworking on "
                          "'eQTL_{}'".format(worker_id, eqtl_index),
                          flush=True)

                # Check if the covariate data still exists. If not, reload it.
                # This hasn't happend yet but I just wanna make sure that if the
                # data is lost from the RAM, the worker will continue.
                try:
                    _ = tech_cov_df.shape
                    _ = cov_df.shape
                except:
                    print("[worker {:2d}]\tcovariate table was lost. "
                          "Reloading data..".format(worker_id), flush=True)
                    tech_cov_df, cov_df = load_covariate_file(cov_inpath,
                                                              tech_covs)

                # Find which order ids to perform.
                order_ids_to_perform = my_work[eqtl_index]

                # Get the complete genotype row for the permutation later.
                genotype_all = geno_df.iloc[row_index, :].copy()

                # Get the missing genotype indices.
                indices = np.arange(geno_df.shape[1])
                eqtl_indices = indices[~geno_df.iloc[row_index, :].isnull().values]

                # Subset the row and present samples for this eQTL.
                genotype = geno_df.iloc[row_index, eqtl_indices].copy()
                expression = expr_df.iloc[row_index, eqtl_indices].copy()
                technical_covs = tech_cov_df.iloc[:, eqtl_indices].copy()
                covariates = cov_df.iloc[:, eqtl_indices].copy()

                # Check if SNP index are identical.
                if genotype.name != expression.name:
                    print("[worker {:2d}]\tindices do "
                          "not match".format(worker_id), flush=True)
                    continue

                # Create the null model. Null model are all the technical
                # covariates multiplied with the genotype + the SNP.
                tech_inter_matrix = technical_covs.mul(genotype, axis=1)
                tech_inter_matrix.index = ["{}_X_SNP".format(x) for x in technical_covs.index]
                intercept = pd.DataFrame(1, index=genotype.index, columns=["intercept"])
                base_matrix = reduce(lambda left, right: pd.merge(left,
                                                                  right,
                                                                  left_index=True,
                                                                  right_index=True),
                                     [intercept,
                                      genotype.to_frame(),
                                      technical_covs.T,
                                      tech_inter_matrix.T])

                # Initialize results variable.
                results = {x: [(genotype.name, genotype.name)] for x in order_ids_to_perform}

                # Loop over the covariates.
                for cov_index in range(cov_df.shape[0]):
                    # Get the covariate we are processing.
                    covarate = covariates.iloc[cov_index, :]
                    cov_name = covarate.name

                    # if verbose:
                    #     print("[worker {:2d}]\tworking on "
                    #           "'{}'".format(worker_id, cov_name),
                    #           flush=True)

                    # Add the covariate to the null matrix if it isn't already.
                    null_matrix = base_matrix.copy()
                    if cov_name not in null_matrix.columns:
                        null_matrix = null_matrix.merge(covarate.to_frame(),
                                                        left_index=True,
                                                        right_index=True)

                    # Create the null model.
                    n_null = null_matrix.shape[0]
                    df_null, rss_null, _ = create_model(null_matrix, expression)

                    # Loop over each permutation sample order. The first order
                    # is the normal order and the remainder are random shuffles.
                    for order_id, sample_order in enumerate(sample_orders):
                        # Making sure I only do the order that I have to do. It
                        # might be the case I got assigned an incomplete eQTL
                        # (i.e. not all sample orders have to be performed) due
                        # to another core dieing halfway through the analysis.
                        if order_id not in order_ids_to_perform:
                            continue

                        # if verbose:
                        #     print("[worker {:2d}]\tworking on "
                        #           "'order_{}'".format(worker_id, order_id),
                        #           flush=True)

                        # Send an heartbeat back to the manager.
                        now = int(time.time())
                        if (now - last_message) > update_frequency:
                            if verbose:
                                print("[worker {:2d}]\tSending "
                                      "heartbeat".format(worker_id), flush=True)
                            result_q.put(
                                (worker_id, "working", None, None, None))
                            last_message = now

                        # Reorder the covariate based on the sample order.
                        # Make sure the labels are in the same order, just
                        # shuffle the values.
                        covariate_all = cov_df.iloc[cov_index, :].copy()
                        covariate_all_index = covariate_all.index
                        covariate_all = covariate_all.reindex(covariate_all.index[sample_order])
                        covariate_all.index = covariate_all_index

                        # Calculate the interaction effect of the covariate of
                        # interest. Then drop the na's from the interaction
                        # term.
                        inter_of_interest = covariate_all * genotype_all
                        inter_name = "{}_X_SNP".format(cov_name)
                        if inter_name in null_matrix.columns:
                            inter_name = inter_name + "_2"
                        inter_of_interest.name = inter_name
                        inter_of_interest = inter_of_interest.iloc[eqtl_indices]

                        # Check if the drop is identical (see above).
                        if not inter_of_interest.index.equals(null_matrix.index):
                            print("[worker {:2d}]\terror in permutation "
                                  "reordering (ID: {})".format(worker_id,
                                                               order_id),
                                  flush=True)
                            results[order_id].append((np.nan, np.nan))
                            continue

                        # Create the alternative matrix and add the interaction
                        # term.
                        alt_matrix = null_matrix.copy()
                        alt_matrix = alt_matrix.merge(inter_of_interest.to_frame(),
                                                      left_index=True,
                                                      right_index=True)

                        # Create the alternative model.
                        n_alt = alt_matrix.shape[0]
                        df_alt, rss_alt, tvalue = create_model(alt_matrix, expression, inter_name)

                        # Make sure the n's are identical.
                        if n_null != n_alt:
                            print("[worker {:2d}]\terror due to unequal "
                                  "n_null and n_alt".format(worker_id),
                                  flush=True)
                            results[order_id].append((np.nan, np.nan))
                            continue

                        # Compare the null and alternative model.
                        pvalue = compare_models(rss_null, rss_alt,
                                                df_null, df_alt,
                                                n_null)

                        # Safe the pvalue.
                        results[order_id].append((pvalue, tvalue))

                # Push the result back to the manager.
                count = 0
                for oid, data in results.items():
                    if len(data) == (cov_df.shape[0] + 1):
                        count += 1
                        result_q.put((worker_id, "result", eqtl_index, oid, data))
                        time.sleep(0.1)
                if verbose:
                    print("[worker {:2d}]\tpushed {} "
                          "results".format(worker_id, count))

            # Wait before getting a new job.
            time.sleep(sleep_time * 2.1)

        except queue.Empty:
            break


def load_covariate_file(cov_inpath, tech_covs):
    """
    Method for loading the covariate file and splitting it into technical
    covariates and covariates of intrest.

    :param cov_inpath: string, the input path of the covariate input file.
    :param tech_covs: list, the indices of the technical covariates.
    :return tech_cov_df: DataFrame, the technical covariates.
    :return cov_df: DataFrame, the covariates of interest.
    """
    cov_df = pd.read_csv(cov_inpath,
                         sep="\t",
                         header=0,
                         index_col=0)

    tech_cov_df = cov_df.loc[tech_covs, :]

    return tech_cov_df, cov_df


def create_model(X, y, name=None):
    """
    Method for creating a multilinear model.

    :param X: DataFrame, the matrix with rows as samples and columns as
                         dimensions.
    :param y: Series, the outcome values.
    :param name: string, the name of the variable accessing for the t-value.
    :return df: int, the degrees of freedom of this model.
    :return ssr: float, the residual sum of squares of this fit.
    :return tvalue: float, beta  / std error.
    """
    # Filter on features with no standard deviation.
    stds = X.std(axis=0) > 0
    stds['intercept'] = True
    X = X.loc[:, stds]

    # Perform the Ordinary least squares fit.
    ols = sm.OLS(y.values, X)
    ols_result = ols.fit()

    df = ols_result.df_model
    ssr = ols_result.ssr

    tvalue = 0
    if name and name in X.columns:
        coef = ols_result.params[name]
        std_err = ols_result.bse[name]
        if std_err > 0:
            tvalue = coef / std_err

    return df, ssr, tvalue


def compare_models(rss1, rss2, df1, df2, n):
    """
    Method for comparing the two models using an F-distribution.

    :param rss1: float, the residual sum of squares of the null model.
    :param rss2: float, the residual sum of squares of the alternative model.
    :param df1: int, the degrees of freedom of the null model.
    :param df2: int, the degrees of freedom of the alternative model.
    :param n: int, the number of samples in the model.
    :return p_value: float, the p-value of the comparison.
    """
    f_value = calc_f_value(rss1, rss2, df1, df2, n)
    p_value = get_p_value(f_value, df1, df2, n)

    return p_value


def calc_f_value(rss1, rss2, df1, df2, n):
    """
    Method for comparing the risdual sum squared of two models using
    the F statistic.

    :param rss1: float, the residual sum of squares of the null model.
    :param rss2: float, the residual sum of squares of the alternative model.
    :param df1: int, the degrees of freedom of the null model.
    :param df2: int, the degrees of freedom of the alternative model.
    :param n: int, the number of samples in the model.
    :return : float, the p-value of the comparison.
    """
    if df1 >= df2:
        return np.nan
    if rss2 >= rss1:
        return 0
    if df2 >= n:
        return 0

    return ((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))


def get_p_value(f_value, df1, df2, n):
    """
    Method for getting the p-value corresponding to a F-distribution.

    :param f_value: float, the f-value.
    :param df1: int, the degrees of freedom of the null model.
    :param df2: int, the degrees of freedom of the alternative model.
    :param n: int, the number of samples in the model.
    :return : float, the p-value corresponding to the f-value.
    """
    # Lower and upper limit of stats.f.sf
    # stats.f.sf(1827.95, dfn=1, dfd=3661) = 5e-324
    # stats.f.sf(9.9e-12, dfn=1, dfd=3661) = 0.9999974567714613

    # Lower and upper limit of cdf
    # stats.f.cdf(69, dfn=1, dfd=3661) = 0.9999999999999999
    # stats.f.cdf(1e-320, dfn=1, dfd=3661) = 1.0730071046473278e-160

    return stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2))
