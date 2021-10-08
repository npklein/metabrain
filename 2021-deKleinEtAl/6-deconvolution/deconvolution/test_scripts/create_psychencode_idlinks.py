#!/usr/bin/env python3

"""
File:         create_psychencode_idlinks.py
Created:      2021/10/08
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
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Create PsychENCODE ID Links"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@st.hanze.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


"""
Syntax:
./create_psychencode_idlinks.py
"""


class main():
    def __init__(self):
        self.phenotype_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-03-09.brain.phenotypes.txt"
        self.psychencode_cf_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/DER-24_Cell_fractions_Normalized.xlsx"
        self.std_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexEUR-cis-Normalised/data/SampleToDataset.txt.gz"

    def start(self):
        print("Step 1: loading data")
        phenotype_df = self.load_file(path=self.phenotype_path, header=0, index_col=None, low_memory=False)
        psychencode_cf_df = self.load_file(path=self.psychencode_cf_path, sheet_name="Sheet1")
        std_df = self.load_file(path=self.std_path, header=0, index_col=None, low_memory=False)
        metabrain_samples = set(std_df["sample"])

        # Pre-process.
        phenotype_df = phenotype_df.loc[:, ["rnaseq_id", "genotype_id", "SampleFull"]]
        phenotype_df.fillna("", inplace=True)
        phenotype_df["combined_id"] = (phenotype_df["rnaseq_id"] + "|" + phenotype_df["genotype_id"] + "|" + phenotype_df["SampleFull"]).str.upper()
        combined_ids = set(phenotype_df["combined_id"])
        phenotype_df.set_index("combined_id", inplace=True)
        print(phenotype_df)

        # Check if substrings in psychencode ids.
        psychencode_ids = psychencode_cf_df.columns.tolist()
        sample_substrings = {}
        for sample1 in psychencode_ids:
            substrings = []
            for sample2 in psychencode_ids:
                if (sample2 != sample1) and ((sample2 in sample1) or (sample1 in sample2)):
                    substrings.append(sample2)

            if len(substrings) > 0:
                sample_substrings[sample1] = substrings

        print("Step 3: Matching")
        psychencode_matches = []
        found_rnaseq_ids = set()
        doubles = []
        n_matched = 0
        for psychencode_id in psychencode_cf_df.columns:
            psychencode_id_str_upper = str(psychencode_id).upper()

            match_info = [psychencode_id, None, None, None]
            for combined_id in combined_ids:
                if psychencode_id_str_upper in combined_id:
                    # Check for substrings.
                    skip = False
                    if psychencode_id in sample_substrings:
                        for sample_substring in sample_substrings[psychencode_id]:
                            if sample_substring in combined_id:
                                skip = True
                                break

                    if skip:
                        continue

                    rnaseq_id, genotype_id, sample_full = phenotype_df.loc[combined_id, :]
                    rnaseq_id_str_upper = str(rnaseq_id).upper()
                    genotype_id_str_upper = str(genotype_id).upper()
                    sample_full_str_upper = str(sample_full).upper()

                    found = False
                    if psychencode_id_str_upper == rnaseq_id_str_upper:
                        match_info[2] = "rnaseq_id"
                        match_info[3] = rnaseq_id
                        found = True
                    elif psychencode_id_str_upper == genotype_id_str_upper:
                        match_info[2] = "genotype_id"
                        match_info[3] = genotype_id
                        found = True
                    elif psychencode_id_str_upper in sample_full_str_upper:
                        match_info[2] = "sample_full"
                        match_info[3] = sample_full
                        found = True
                    else:
                        pass

                    if found:
                        match_info[1] = rnaseq_id
                        n_matched += 1
                        break

            if match_info[1] in found_rnaseq_ids:
                doubles.append(match_info[1])

            if match_info[1] is not None:
                found_rnaseq_ids.add(match_info[1])

            psychencode_matches.append(match_info)

        link_df = pd.DataFrame(psychencode_matches, columns=["PsychENCODE", "MetaBrain", "ID name", "ID value"])
        print(link_df)
        print(link_df["ID name"].value_counts())
        print("\tN-matched {}/{} [{:.2f}%]".format(n_matched, psychencode_cf_df.shape[1], (100 / psychencode_cf_df.shape[1]) * n_matched))

        link_df.dropna(inplace=True)
        self.save_file(df=link_df, outpath="PsychENCODE_ID_links.txt.gz", index=False)

        for double in doubles:
            print("Double: {}".format(double))
            print(link_df.loc[link_df["MetaBrain"] == double, :])
            print("")

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None,
                  skiprows=0, sheet_name=None, low_memory=True):
        if path.endswith(".xlsx"):
            df = pd.read_excel(path, header=header, index_col=index_col,
                         nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        else:
            df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                             nrows=nrows, skiprows=skiprows,
                             low_memory=low_memory)

        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))


if __name__ == '__main__':
    m = main()
    m.start()
