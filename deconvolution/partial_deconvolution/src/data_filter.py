"""
File:         data_filter.py
Created:      2020/09/04
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

# Third party imports.
import pandas as pd

# Local application imports.


class DataFilter:
    def __init__(self, settings, raw_expression):
        self.sample_annotation_path = settings.get_sample_annotation_path()
        self.sample_id = settings.get_sample_id()
        self.sample_filter_path = settings.get_sample_filter_path()
        self.cohort_id = settings.get_cohort_id()
        self.cohort_filter = settings.get_cohort_filter()
        self.annotation_id = settings.get_annotation_id()
        self.annotation_filter = settings.get_annotation_filter()
        self.cohort_corr = settings.get_cohort_corr()
        self.raw_expression = raw_expression

        self.sample_annotation = None
        self.expression = None
        self.cohorts = None

    def work(self):
        # Load the sample annotation file
        self.sample_annotation = self.load_sample_annotation()

        # Filter the expression.
        self.expression = self.filter_expression_data(self.raw_expression)

        # Create the cohorts matrix.
        if self.cohort_corr:
            self.cohorts = self.create_cohorts_matrix()

    def load_sample_annotation(self):
        print("Loading samples annotation")
        return pd.read_csv(self.sample_annotation_path,
                           sep="\t",
                           low_memory=False)

    def get_translate_dict(self, value, key=None):
        if key is None:
            key = self.sample_id
        if self.sample_annotation is None or value is None:
            return None

        return dict(zip(self.sample_annotation.loc[:, key],
                        self.sample_annotation.loc[:, value]))

    def filter_expression_data(self, raw_expression):
        print("Filtering expression data")

        sample_filter = self.load_sample_filter()
        cohort_dict = self.get_translate_dict(key=self.sample_id,
                                              value=self.cohort_id)
        annotation_dict = self.get_translate_dict(key=self.sample_id,
                                                  value=self.annotation_id)

        include = []
        for sample in raw_expression.columns:
            if (sample in sample_filter) and \
                (cohort_dict is None or sample in cohort_dict.keys()) and \
                (cohort_dict is None or self.cohort_filter is None or str(cohort_dict[sample]).upper() in self.cohort_filter) and \
                    (annotation_dict is None or self.annotation_filter is None or str(annotation_dict[sample]).upper() in self.annotation_filter):
                include.append(sample)

        return raw_expression.loc[:, include]

    def load_sample_filter(self):
        print("Loading sample filter")
        if self.sample_filter_path is None:
            return None
        else:
            df = pd.read_csv(self.sample_filter_path, sep="\t",
                             header=0, index_col=0,
                             nrows=1)
            return list(df.columns)

    def create_cohorts_matrix(self):
        print("Creating cohort matrix.")

        sample_to_cohort_dict = self.get_translate_dict(key=self.sample_id,
                                                        value=self.cohort_id)

        cohort_to_sample_dict = {}
        for sample in self.expression.columns:
            cohort = sample_to_cohort_dict[sample]
            if cohort in cohort_to_sample_dict.keys():
                value = cohort_to_sample_dict[cohort]
                value.append(sample)
                cohort_to_sample_dict[cohort] = value
            else:
                cohort_to_sample_dict[cohort] = [sample]

        cohort_df = pd.DataFrame(0,
                                 index=sample_to_cohort_dict.keys(),
                                 columns=cohort_to_sample_dict.keys())
        for cohort in cohort_df.columns:
            cohort_df.loc[cohort_df.index.isin(cohort_to_sample_dict[cohort]), cohort] = 1

        return cohort_df

    def get_expression(self):
        return self.expression

    def get_cohorts(self):
        return self.cohorts

    def get_shape_diff(self):
        return self.raw_expression.shape[0] - self.expression.shape[0], self.raw_expression.shape[1] - self.expression.shape[1]

    def print_info(self):
        print("Pre-filter shape:")
        print("\t{} rows and {} columns".format(self.raw_expression.shape[0], self.raw_expression.shape[1]))
        print("Post-filter shape:")
        print("\t{} rows and {} columns".format(self.expression.shape[0], self.expression.shape[1]))
        print("Shape difference:")
        print("\t{} rows and {} columns".format(self.get_shape_diff()[0], self.get_shape_diff()[1]))
        print("")
