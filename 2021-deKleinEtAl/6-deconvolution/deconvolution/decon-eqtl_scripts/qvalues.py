"""
File:         qvalues.py
Created:      2020/07/21
Last Changed:
Author:       M.Vochteloo

Python implementation of R Bioconductor’s qvalue package Version 2.24.0.
"""

# Standard imports.

# Third party imports.
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from scipy.stats import norm


# Local application imports.


def qvalue(self, p, fdr_level=None, pfdr=False, lfdr_out=True, pi0=None,
           **kwargs):
    p_in = np.copy(p)
    qvals_out = np.copy(p)
    lfdr_out_a = np.copy(p)
    rm_na = ~np.isnan(p)
    p = p[rm_na]
    if np.min(p) < 0 or np.max(p) > 1:
        print("p-values not in valid range [0, 1].")
        return None
    if fdr_level is not None and (fdr_level <= 0 or fdr_level > 1):
        print("'fdr.level' must be in (0, 1].")
        return None
    if pi0 is not None and (pi0 <= 0 or pi0 > 1):
        print("pi0 is not (0,1]")
        return None

    # Estimate pi0
    pi0s = {"pi0": pi0,
            "pi0_lambda": None,
            "lambda_values": None,
            "pi0_smooth": None}
    if pi0 is None:
        pi0s = self.pi0est(p=p, **kwargs)
        pi0 = pi0s["pi0"]

    m = np.size(p)
    i = np.arange(m, 0, -1)
    o = np.argsort(p)[::-1]
    ro = np.argsort(o)

    if pfdr:
        rank = i * (1 - (1 - p[o]) ^ m)
    else:
        rank = i
    qvals = pi0 * np.min(
        (np.ones_like(p), pd.Series(p[o] * m / rank).cummin().to_numpy()),
        axis=0)[ro]
    qvals_out[rm_na] = qvals

    if lfdr_out:
        lfdr_out_a[rm_na] = lfdr(p=p, pi0=pi0, **kwargs)
    else:
        lfdr_out_a = None

    retval = {
        "pi0": pi0,
        "qvalues": qvals_out,
        "pvalues": p_in,
        "lfdr_out": lfdr_out_a,
        "pi0.lambda": pi0s["pi0_lambda"],
        "lambda": pi0s["lambda_values"],
        "pi0.smooth": pi0s["pi0_smooth"]
    }
    if fdr_level is not None:
        retval["fdr_level"] = fdr_level
        retval["significant"] = qvals <= fdr_level

    return retval


def lfdr(self, p, pi0, trunc=True, monotone=True, transf="logit", adj=1.5,
         eps=1e-08):
    if transf not in ["probit", "logit"]:
        print("ERROR: transf invalid")
        return None
    lfdr_out = np.copy(p)
    rm_na = ~np.isnan(p)
    p = p[rm_na]

    if np.min(p) < 0 or np.max(p) > 1:
        print("P-values not in valid range [0,1].")
        return lfdr_out

    eps_a = np.full_like(p, eps)
    if transf == "probit":
        p = np.max((p, eps_a), axis=0)
        p = np.min((p, 1 - eps_a), axis=0)
        x = norm.ppf(p)
        mydx, mydy = self.density(x, adjust=adj)
        y = self.smooth_spline(x_train=mydx, y_train=mydy, x_test=x)
        lfdr_values = pi0 * norm.pdf(x) / y
    elif transf == "logit":
        x = np.log((p + eps_a) / (1 - p + eps_a))
        mydx, mydy = self.density(x, adjust=adj)
        y = self.smooth_spline(x_train=mydx, y_train=mydy, x_test=x)
        dx = np.exp(x) / np.power((1 + np.exp(x)), 2)
        lfdr_values = (pi0 * dx) / y
    else:
        print("ERROR: transf must be one of \"probit\" or \"logit\".")
        return lfdr_out

    if trunc:
        lfdr_values[lfdr_values > 1] = 1

    if monotone:
        o = np.argsort(p)[::-1]
        ro = np.argsort(o)
        lfdr_values = pd.Series(lfdr_values[o]).cummax().to_numpy()[ro]

    lfdr_out[rm_na] = lfdr_values

    return lfdr_out


def pi0est(self, p, lambda_values=np.arange(0.05, 1, 0.05),
           pio_method="smoother",
           smooth_df=3, smooth_log_pi0=False):
    retval = {"pi0": None,
              "pi0_lambda": None,
              "lambda_values": None,
              "pi0_smooth": None}
    if pio_method not in ["smoother", "bootstrap"]:
        print("ERROR: pio method invalid")
        return retval
    if not type(lambda_values) is np.ndarray:
        lambda_values = np.array([lambda_values])

    rm_na = ~np.isnan(p)
    p = p[rm_na]
    m = np.size(p)
    lambda_values = np.sort(lambda_values)
    ll = np.size(lambda_values)

    retval["lambda_values"] = lambda_values

    if np.min(p) < 0 or np.max(p) > 1:
        print("ERROR: p-values not in valid range [0, 1].")
        return retval
    if 1 < ll < 4:
        print("ERROR: length(lambda)= {}. If length of lambda greater than"
              " 1, you need at least 4 values.".format(ll))
        return retval
    if np.min(lambda_values) < 0 or np.max(lambda_values) > 1:
        print("ERROR: Lambda must be within [0, 1).")
        return retval

    if ll == 1:
        pi0 = np.mean(p >= lambda_values[0]) / (1 - lambda_values[0])
        retval["pi0_lambda"] = pi0
        pi0 = min(pi0, 1)
    else:
        bins = np.digitize(p, lambda_values)
        bincounts = dict(zip(*np.unique(bins, return_counts=True)))
        frequencies = np.array([bincounts[x] if x in bincounts else 0 for x in
                                np.arange(ll, 0, -1)])
        pi0 = np.cumsum(frequencies) / (m * (1 - lambda_values[::-1]))
        pi0 = pi0[::-1]
        retval["pi0_lambda"] = pi0
        if pio_method == "smoother":
            if smooth_log_pi0:
                pi0 = np.log(pi0)
            pi0_smooth = np.exp(self.smooth_spline(x_train=lambda_values,
                                                   y_train=pi0,
                                                   df=smooth_df))
            pi0 = min(pi0_smooth[ll - 1], 1)
            retval["pi0_smooth"] = pi0_smooth
        elif pio_method == "bootstrap":
            minpi0 = np.quantile(pi0, q=0.1)
            W = np.array([np.sum(p >= l) for l in lambda_values],
                         dtype=np.float64)
            mse = (W / (m * m * np.power((1 - lambda_values), 2))) * (
                        1 - W / m) + np.power((pi0 - minpi0), 2)
            pi0 = min(pi0[mse == np.min(mse)][0], 1)
        else:
            print("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
            return retval

    if pi0 <= 0:
        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")

    retval["pi0"] = pi0

    return retval


def smooth_spline(x_train, y_train, x_test=None, **kwargs):
    r_x_train = robjects.FloatVector(x_train)
    r_y_train = robjects.FloatVector(y_train)
    r_x_test = r_x_train
    if x_test is not None:
        r_x_test = robjects.FloatVector(x_test)
    spline_xy = robjects.r['smooth.spline'](r_x_train, r_y_train, **kwargs)
    return np.array(robjects.r['predict'](spline_xy, r_x_test).rx2('y'))


def density(x, adjust):
    r_x = robjects.FloatVector(x)
    density_func = robjects.r['density'](r_x, adjust=adjust)
    return np.array(density_func.rx2('x')), np.array(density_func.rx2('y'))