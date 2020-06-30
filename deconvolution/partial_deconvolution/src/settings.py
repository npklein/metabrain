"""
File:         settings.py
Created:      2020/06/29
Last Changed: 2020/06/30
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
from datetime import datetime
import json
import os

# Third party imports.

# Local application imports.


class Settings:
    def __init__(self, data_path, signature_path, translate_path, sample_path,
                 cohort, ground_truth_path, min_expr, normalize, zscore, log2,
                 decon_method, sum_to_one, extension):
        self.data_path = data_path
        self.signature_path = signature_path
        self.translate_path = translate_path
        self.sample_path = sample_path
        self.cohort = cohort
        self.ground_truth_path = ground_truth_path
        self.min_expr = min_expr
        self.normalize = normalize
        self.zscore = zscore
        self.log2 = log2
        self.decon_method = decon_method
        self.sum_to_one = sum_to_one
        self.extension = extension

        self.outpath = None
        self.real_info_per_celltype = None
        self.filter_shape_diff = None
        self.sign_shift = None
        self.expr_shift = None
        self.avg_residuals = None
        self.pred_info_per_celltype = None
        self.comparison_n_samples = None
        self.comparison_rss = None

    def get_data_path(self):
        return self.data_path

    def get_signature_path(self):
        return self.signature_path

    def get_translate_path(self):
        return self.translate_path

    def get_sample_path(self):
        return self.sample_path

    def get_cohort(self):
        return self.cohort

    def get_ground_truth_path(self):
        return self.ground_truth_path

    def get_ground_truth_type(self):
        return os.path.basename(self.ground_truth_path).split(".")[0].replace("_counts", "").replace("_", " ")

    def get_min_expr(self):
        return self.min_expr

    def get_normalize(self):
        return self.normalize

    def get_zscore(self):
        return self.zscore

    def get_log2(self):
        return self.log2

    def get_decon_method(self):
        return self.decon_method

    def get_sum_to_one(self):
        return self.sum_to_one

    def get_extension(self):
        return self.extension

    def set_output_path(self, outpath):
        self.outpath = outpath

    def get_output_path(self):
        return self.outpath

    def set_real_info_per_celltype(self, real_info_per_celltype):
        self.real_info_per_celltype = real_info_per_celltype

    def set_filter_shape_diff(self, filter_shape_diff):
        self.filter_shape_diff = filter_shape_diff

    def set_sign_shift(self, sign_shift):
        self.sign_shift = sign_shift

    def set_expr_shift(self, expr_shift):
        self.expr_shift = expr_shift

    def set_avg_residuals(self, avg_residuals):
        self.avg_residuals = avg_residuals

    def set_pred_info_per_celltype(self, pred_info_per_celltype):
        self.pred_info_per_celltype = pred_info_per_celltype

    def set_comparison_n_samples(self, comparison_n_samples):
        self.comparison_n_samples = comparison_n_samples

    def set_comparison_rss(self, comparison_rss):
        self.comparison_rss = comparison_rss

    def get_title(self):
        return "{} partial deconvolution using {}".format(self.cohort, self.decon_method)

    def get_subtitle(self):
        avg_resid_str = ""
        if self.avg_residuals is not None:
            avg_resid_str = "{:.2f}".format(self.avg_residuals)

        rss_str = ""
        if self.comparison_rss is not None:
            rss_str = "{:.2f}".format(self.comparison_rss)

        return "min.expr.: {}, norm.: {}, z-score: {}, log2: {}, " \
               "sum-to-one: {}\n avg.residuals: {}, N: {}, RSS: {}".format(
            self.min_expr, self.normalize, self.zscore, self.log2,
            self.sum_to_one, avg_resid_str, self.comparison_n_samples, rss_str)

    def save_settings(self):
        if self.outpath is None:
            print("No output path is set.")
            return

        data = {"date time": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                "data_path": self.data_path,
                "signature_path": self.signature_path,
                "translate_path": self.translate_path,
                "sample_path": self.sample_path,
                "cohort": self.cohort,
                "ground_truth_path": self.ground_truth_path,
                "min_expr": self.min_expr,
                "normalize": self.normalize,
                "zscore": self.zscore,
                "log2": self.log2,
                "decon_method": self.decon_method,
                "sum_to_one": self.sum_to_one,
                "extension": self.extension,
                "real_info_per_celltype": self.real_info_per_celltype,
                "filter_shape_diff": self.filter_shape_diff,
                "sign_shift": self.sign_shift,
                "expr_shift": self.expr_shift,
                "avg_residuals": self.avg_residuals,
                "pred_info_per_celltype": self.pred_info_per_celltype,
                "comparison_n_samples": self.comparison_n_samples,
                "comparison_rss": self.comparison_rss}

        outpath = os.path.join(self.outpath, 'settings.json')
        with open(outpath, 'w') as f:
            json.dump(data, f)
        f.close()
