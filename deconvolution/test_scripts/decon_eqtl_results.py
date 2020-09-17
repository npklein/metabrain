#!/usr/bin/env python3

"""
File:         decon_eqtl_results.py
Created:      2020/09/16
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
import os

# Third party imports.
import pandas as pd
from scipy import stats
from statsmodels.stats import multitest

# Local application imports.

# Metadata
__program__ = "Decon eQTL Results"
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


class main():
    def __init__(self):
        self.decon_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/2020-07-16-decon-eQTL/cis/cortex/decon_out/deconvolutionResults.csv"
        self.gene_info_path = "/groups/umcg-biogen/tmp03/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz"
        self.suppl_table_x = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/Supplementary_Table_X_-_LD_overlap.xlsx"
        self.outdir = str(Path(__file__).parent.parent)

    def start(self):
        print("Loading gene info matrix.")
        gene_info_df = pd.read_csv(self.gene_info_path, sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.gene_info_path),
                                      gene_info_df.shape))
        gene_dict = dict(zip(gene_info_df["ArrayAddress"], gene_info_df["Symbol"]))

        print("Loading decon-eQTL cell type interaction results.")
        decon_df = pd.read_csv(self.decon_path, sep="\t", header=0, index_col=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.decon_path),
                                      decon_df.shape))
        print(decon_df)
        decon_df.reset_index(drop=False, inplace=True)
        decon_df = decon_df.melt(id_vars="index", value_vars=[x for x in decon_df.columns if x.endswith("pvalue")])
        decon_df["FDR"] = multitest.multipletests(decon_df['value'], method='fdr_bh')[1]
        decon_df[['ProbeName', 'SNPName']] = decon_df["index"].str.split("_", n=1, expand=True)
        decon_df['HGNCName'] = decon_df["ProbeName"].map(gene_dict)
        decon_df.reset_index(drop=True, inplace=True)
        decon_df.sort_values(by="FDR", ascending=True, inplace=True)
        decon_df = decon_df.loc[decon_df["FDR"] < 0.05, :]
        print(decon_df)

        decon_df.to_csv(os.path.join(self.outdir, "decon-eQTL_cis_signifDeconvolutionResults.txt.gz"), sep="\t", header=True, index=True, compression="gzip")

        print("Loading suppl tabel X.")
        supll_table_df = pd.read_excel(self.suppl_table_x, sheet_name="Cortex-EUR", sep="\t", header=0)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(self.suppl_table_x),
                                      supll_table_df.shape))
        print(supll_table_df)

        df = supll_table_df.merge(decon_df, left_on=["LinkedEQTLSNP", "LinkedEQTLGenes"], right_on=["SNPName", "ProbeName"])
        df.sort_values(by="FDR", ascending=True, inplace=True)
        df.to_csv(os.path.join(self.outdir, "disease_SNP_CTInteractions.txt.gz"), sep="\t", header=True, index=True, compression="gzip")
        print(df)

        cell_types = list(df["variable"].unique())
        interest = ["ieu-a-1239", "2019-MSGWAS", "ebi-a-GCST006572", "ebi-a-GCST006250", "ieu-a-22", "ebi-a-GCST006940"]
        data = []
        indices = []
        for gwas_id in list(df["GWASID"].unique()):
            counts = df.loc[df["GWASID"].str.lower() == gwas_id.lower(), "variable"].value_counts()

            values = []
            for ct in cell_types:
                if ct in counts:
                    values.append(counts[ct])
                else:
                    values.append(0)
            if len(counts.index) > 0 and gwas_id in interest:
                data.append(values)
                indices.append(gwas_id)
        output = pd.DataFrame(data, index=indices, columns=[x.split("_")[0] for x in cell_types])
        print(output)


if __name__ == '__main__':
    m = main()
    m.start()
