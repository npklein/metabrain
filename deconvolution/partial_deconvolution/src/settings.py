"""
File:         settings.py
Created:      2020/06/29
Last Changed: 2020/09/04
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
    def __init__(self, data_path, signature_path, translate_path,
                 ground_truth_path, sample_annotation_path, sample_id,
                 sample_filter_path, cohort_id, cohort_filter, annotation_id,
                 annotation_filter, min_expr, cohort_corr, normalize, zscore,
                 log2, decon_method, sum_to_one, extension):
        self.data_path = data_path
        self.signature_path = signature_path
        self.translate_path = translate_path
        self.ground_truth_path = ground_truth_path
        self.sample_annotation_path = sample_annotation_path
        self.sample_id = sample_id
        self.sample_filter_path = sample_filter_path
        self.cohort_id = cohort_id
        if cohort_filter is not None:
            cohort_filter = [x.replace("_", " ").upper() for x in cohort_filter]
        self.cohort_filter = cohort_filter
        self.annotation_id = annotation_id
        if annotation_filter is not None:
            annotation_filter = [x.replace("_", " ").upper() for x in annotation_filter]
        self.annotation_filter = annotation_filter
        self.min_expr = min_expr
        self.cohort_corr = cohort_corr
        self.normalize = normalize
        self.zscore = zscore
        self.log2 = log2
        self.decon_method = decon_method
        self.sum_to_one = sum_to_one
        self.extension = extension

        self.outdir_path = None
        self.outsubdir_path = None
        self.real_info_per_celltype = None
        self.filter1_shape_diff = None
        self.filter2_shape_diff = None
        self.sign_shift = None
        self.expr_shift = None
        self.n_samples = None
        self.n_genes = None
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

    def get_ground_truth_path(self):
        return self.ground_truth_path

    def get_ground_truth_type(self):
        return os.path.basename(self.ground_truth_path).split(".")[0].replace("_counts", "")

    def get_sample_annotation_path(self):
        return self.sample_annotation_path

    def get_sample_id(self):
        return self.sample_id

    def get_sample_filter_path(self):
        return self.sample_filter_path

    def get_cohort_id(self):
        return self.cohort_id

    def get_cohort_filter(self):
        return self.cohort_filter

    def get_annotation_id(self):
        return self.annotation_id

    def get_annotation_filter(self):
        return self.annotation_filter

    def get_min_expr(self):
        return self.min_expr

    def get_cohort_corr(self):
        return self.cohort_corr

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

    def get_outdir_path(self):
        return self.outdir_path

    def get_outsubdir_path(self):
        return self.outsubdir_path

    def set_outdir_path(self, outdir_path):
        self.outdir_path = outdir_path

    def set_outsubdir_path(self, outsubdir_path):
        self.outsubdir_path = outsubdir_path

    def set_real_info_per_celltype(self, real_info_per_celltype):
        self.real_info_per_celltype = real_info_per_celltype

    def set_filter1_shape_diff(self, filter1_shape_diff):
        self.filter1_shape_diff = filter1_shape_diff

    def set_filter2_shape_diff(self, filter2_shape_diff):
        self.filter2_shape_diff = filter2_shape_diff

    def set_sign_shift(self, sign_shift):
        self.sign_shift = sign_shift

    def set_expr_shift(self, expr_shift):
        self.expr_shift = expr_shift

    def set_n_samples(self, n_samples):
        self.n_samples = n_samples

    def set_n_genes(self, n_genes):
        self.n_genes = n_genes

    def set_avg_residuals(self, avg_residuals):
        self.avg_residuals = avg_residuals

    def set_pred_info_per_celltype(self, pred_info_per_celltype):
        self.pred_info_per_celltype = pred_info_per_celltype

    def set_comparison_n_samples(self, comparison_n_samples):
        self.comparison_n_samples = comparison_n_samples

    def set_comparison_rss(self, comparison_rss):
        self.comparison_rss = comparison_rss

    def get_title(self):
        return "partial deconvolution using {}".format(self.decon_method)

    def get_subtitle(self):
        avg_resid_str = ""
        if self.avg_residuals is not None:
            avg_resid_str = "{:.2f}".format(self.avg_residuals)

        rss_str = ""
        if self.comparison_rss is not None:
            rss_str = "{:.2f}".format(self.comparison_rss)

        return "min.expr.: {}, coh.corr: {}, norm.: {}, z-score: {}, " \
               "log2: {},\nsum-to-one: {}, avg.residuals: {}, N: {}, " \
               "RSS: {}".format(self.min_expr, self.cohort_corr, self.normalize,
                                self.zscore, self.log2, self.sum_to_one,
                                avg_resid_str, self.comparison_n_samples,
                                rss_str)

    def save_data_settings(self):
        data = {"date time": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                "data_path": self.data_path,
                "signature_path": self.signature_path,
                "translate_path": self.translate_path
                }

        self.save_settings(self.outdir_path, 'data_settings', data)

    def save_all_settings(self):
        data = {"date time": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                "data_path": self.data_path,
                "signature_path": self.signature_path,
                "translate_path": self.translate_path,
                "ground_truth_path": self.ground_truth_path,
                "outdir_path": self.outdir_path,
                "sample_annotation_path": self.sample_annotation_path,
                "sample_id": self.sample_id,
                "sample_filter_path": self.sample_filter_path,
                "cohort_id": self.cohort_id,
                "cohort_filter": self.cohort_filter,
                "annotation_id": self.annotation_id,
                "annotation_filter": self.annotation_filter,
                "filter1_shape_diff": self.filter1_shape_diff,
                "cohort_corr": self.cohort_corr,
                "outsubdir_path": self.outsubdir_path,
                "min_expr": self.min_expr,
                "normalize": self.normalize,
                "zscore": self.zscore,
                "log2": self.log2,
                "decon_method": self.decon_method,
                "sum_to_one": self.sum_to_one,
                "extension": self.extension,
                "real_info_per_celltype": self.real_info_per_celltype,
                "filter2_shape_diff": self.filter2_shape_diff,
                "sign_shift": self.sign_shift,
                "expr_shift": self.expr_shift,
                "n_samples": self.n_samples,
                "n_genes": self.n_genes,
                "avg_residuals": self.avg_residuals,
                "pred_info_per_celltype": self.pred_info_per_celltype,
                "comparison_n_samples": self.comparison_n_samples,
                "comparison_rss": self.comparison_rss}
        self.save_settings(self.outsubdir_path, 'all_settings', data)

    @staticmethod
    def save_settings(path, name, data):
        if path is None:
            print("No output path is set.")
            return

        outpath = os.path.join(path, '{}.json'.format(name))
        with open(outpath, 'w') as f:
            json.dump(data, f)
        f.close()
