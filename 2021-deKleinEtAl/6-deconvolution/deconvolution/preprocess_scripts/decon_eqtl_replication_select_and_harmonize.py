#!/usr/bin/env python3

"""
File:         decon_eqtl_replication_select_and_harmonize.py
Created:      2021/06/23
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
import argparse
import gzip
import os

# Third party imports.
import pandas as pd
import numpy as np

# Local application imports.

# Metadata
__program__ = "Decon-eQTL Select and Harmonize"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
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
./decon_eqtl_replication_select_and_harmonize.py -re /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-05-26-eqtls-rsidfix-popfix/cis/2020-05-26-Cortex-AFR/replicateCortex-EUR/eQTLProbesFDR0.05-ProbeLevel.txt.gz -da /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexEUR-cis/create_matrices/genotype_alleles.txt.gz -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_table.txt.gz -al /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/CortexAFR-cis-Replication-EUR/create_matrices/genotype_alleles.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/preprocess_scripts/pre_process_decon_expression_matrix/CortexAFR-cis/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.SampleSelection.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.CovariatesRemovedOLS.ExpAdded.txt.gz -of CortexAFR-cis-Replication-EUR
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.replication_eqtl_path = getattr(arguments, 'replication_eqtl')
        self.discovery_alleles_path = getattr(arguments, 'discovery_alleles')
        self.geno_path = getattr(arguments, 'genotype')
        self.alleles_path = getattr(arguments, 'alleles')
        self.expr_path = getattr(arguments, 'expression')
        outdir = getattr(arguments, 'outdir')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        if outdir is None:
            outdir = str(Path(__file__).parent.parent)
        self.outdir = os.path.join(outdir, 'decon_eqtl_replication_select_and_harmonize', outfolder)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-re",
                            "--replication_eqtl",
                            type=str,
                            required=True,
                            help="The path to the replication eQTLs")
        parser.add_argument("-da",
                            "--discovery_alleles",
                            type=str,
                            required=True,
                            help="The path to the discovery alleles matrix")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-al",
                            "--alleles",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-od",
                            "--outdir",
                            type=str,
                            required=False,
                            default=None,
                            help="The name of the output path.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Loading eQTLs data.")
        eqtl_df = self.load_file(self.replication_eqtl_path, header=0, index_col=None)
        snps_of_interest = eqtl_df["SNPName"]
        genes_of_interest = eqtl_df["ProbeName"]

        print("Loading the allele assessed files.")
        replication_alleles_df = self.load_file(self.alleles_path, header=0, index_col=0)
        discovery_alleles_df = self.load_file(self.discovery_alleles_path, header=0, index_col=0)
        alleles_df = replication_alleles_df[["MinorAllele"]].merge(discovery_alleles_df[["MinorAllele"]], left_index=True, right_index=True)
        alleles_df["flip"] = alleles_df.iloc[:, 0] != alleles_df.iloc[:, 1]
        flip_dict = dict(zip(alleles_df.index, alleles_df["flip"]))
        del replication_alleles_df, discovery_alleles_df

        print("Processing genotype file.")
        geno_df = self.process_file(inpath=self.geno_path,
                                    filter_set=set(snps_of_interest),
                                    flip_dict=flip_dict)

        print("\tReordering matrix.")
        geno_df = geno_df.loc[snps_of_interest, :]
        geno_df.index.name = None
        mask = geno_df.isnull().all(1).to_numpy()
        geno_df = geno_df.loc[~mask, :]

        print("\tSaving output file.")
        self.save_file(df=geno_df,
                       outpath=os.path.join(self.outdir, os.path.basename(self.geno_path).replace(".gz", "")))
        del geno_df

        print("Processing expression file.")
        expr_df = self.process_file(inpath=self.expr_path,
                                    filter_set=set(genes_of_interest))

        print("\tReordering matrix.")
        expr_df = expr_df.loc[genes_of_interest, :]
        expr_df.index.name = None
        expr_df = expr_df.loc[~mask, :]

        print("\tSaving output file.")
        self.save_file(df=expr_df,
                       outpath=os.path.join(self.outdir, os.path.basename(self.expr_path).replace(".gz", "")))
        del expr_df

        print("Processing SnpsToTest file.")
        snps_to_test_df = eqtl_df[["ProbeName", "SNPName"]].copy()
        snps_to_test_df = snps_to_test_df.loc[~mask, :]
        self.save_file(df=snps_to_test_df,
                       outpath=os.path.join(self.outdir, "snp_gene_list.txt"),
                       index=False)

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def process_file(inpath, filter_set, flip_dict=None, sep="\t",
                     header=0, index_col=0):
        found_indices = set()
        info = []
        with gzip.open(inpath, 'rb') as f:
            for i, line in enumerate(f):
                splitted_line = line.decode().strip('\n').split(sep)
                if i == 0:
                    info.append(splitted_line)
                else:
                    index = splitted_line[0]
                    if index not in filter_set or index in found_indices:
                        continue

                    if flip_dict is not None and flip_dict[index]:
                        values = [float(x) if x != '' else np.nan for x in splitted_line[1:]]
                        data = np.array(values)
                        data = 2 - data
                        info.append([index] + data.tolist())
                    else:
                        info.append(splitted_line)

                    found_indices.add(index)
        f.close()

        df = pd.DataFrame(info)
        df.index = df.iloc[:, index_col]
        df.index.name = "-"
        df = df.drop(df.columns[index_col], axis=1)
        df.columns = df.iloc[header, :]
        df = df.drop(df.index[header], axis=0)
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

    def print_arguments(self):
        print("Arguments:")
        print("  > Replication eQTL path: {}".format(self.replication_eqtl_path))
        print("  > Discovery alleles path: {}".format(self.discovery_alleles_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Alleles path: {}".format(self.alleles_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Output directory: {}".format(self.outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
