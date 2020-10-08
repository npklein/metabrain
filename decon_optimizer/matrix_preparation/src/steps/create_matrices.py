"""
File:         create_matrices.py
Created:      2020/10/08
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
import gzip
import os

# Third party imports.
import numpy as np
import pandas as pd

# Local application imports.
from utilities import prepare_output_dir, check_file_exists, load_dataframe, save_dataframe, construct_dict_from_df


class CreateMatrices:
    def __init__(self, settings, sample_dict, sample_order, eqtl_df,
                 force, outdir):
        self.geno_file = settings["genotype_datafile"]
        self.expr_file = settings["expression_datafile"]
        self.sign_file = settings["signature_datafile"]
        gene_translate_settings = settings["gene_translate"]
        self.gene_info_file = gene_translate_settings["datafile"]
        self.ens_id = gene_translate_settings["ENS_column"]
        self.hgnc_id = gene_translate_settings["HGNC_column"]
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.eqtl_df = eqtl_df
        self.force = force
        self.print_interval = 250

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_matrices')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.geno_outpath = os.path.join(self.outdir, "genotype_table.txt.gz")
        self.alleles_outpath = os.path.join(self.outdir, "genotype_alleles.txt.gz")
        self.expr_outpath = os.path.join(self.outdir, "expression_table.txt.gz")
        self.sign_expr_outpath = os.path.join(self.outdir, "signature_expr_table.txt.gz")

        # Create empty variable.
        self.sign_df = None
        self.gene_info_df = None
        self.geno_df = None
        self.alleles_df = None
        self.expr_df = None
        self.sign_expr_df = None

    def start(self):
        print("Starting creating matrices.")
        self.print_arguments()

        if not check_file_exists(self.geno_outpath) or not check_file_exists(self.alleles_outpath) or self.force:
            print("Parsing genotype input data.")
            alleles_df, genotype_df = self.parse_genotype_file()

            print("Filter and save.")
            for df, outpath in zip([genotype_df, alleles_df], [self.geno_outpath, self.alleles_outpath]):
                df = self.fill_and_reorder(df, "SNPName")
                save_dataframe(df=df, outpath=outpath,
                               index=True, header=True)

        if not check_file_exists(self.expr_outpath) or not check_file_exists(self.sign_expr_outpath) or self.force:
            print("Loading signature matrix.")
            self.sign_df = load_dataframe(inpath=self.sign_file,
                                          header=0,
                                          index_col=0)
            signature_genes = set(self.sign_df.index.to_list())

            print("Loading gene traslate dict.")
            self.gene_info_df = load_dataframe(inpath=self.gene_info_file,
                                               header=0,
                                               index_col=None)
            gene_trans_dict = construct_dict_from_df(self.gene_info_df, self.ens_id, self.hgnc_id)

            print("Parsing expression data.")
            expression_df, signature_expr_df = self.parse_expression_file(signature_genes, gene_trans_dict)

            print("Filter and save.")
            for df, outpath in zip([expression_df, signature_expr_df], [self.expr_outpath, self.sign_expr_outpath]):
                df = self.fill_and_reorder(df, "ProbeName")
                save_dataframe(df=df, outpath=outpath,
                               index=True, header=True)

    def parse_genotype_file(self):
        interest = set(self.eqtl_df.loc[:, "SNPName"].to_list())

        alleles_data_collection = []
        genotype_data_collection = []
        indices = []
        with gzip.open(self.geno_file, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    print("\t\tprocessed {} lines\tfound {}/{} genotype lines.".format(i, len(indices), self.eqtl_df.shape[0]))

                splitted_line = np.array(line.decode().strip('\n').split('\t'))
                index = splitted_line[0]
                data = splitted_line[1:]
                if i == 0:
                    columns = data
                else:
                    alleles_data = data[:2]
                    genotype_data = splitted_line[2:]
                    if index in interest:
                        indices.append(index)
                        alleles_data_collection.append(alleles_data)
                        genotype_data_collection.append([float(x) for x in genotype_data])

        f.close()

        alleles_df = pd.DataFrame(alleles_data_collection, index=indices, columns=columns)
        genotype_df = pd.DataFrame(genotype_data_collection, index=indices, columns=columns)

        return alleles_df, genotype_df

    def parse_expression_file(self, signature_genes, gene_trans_dict):
        expression_interest = set(self.eqtl_df.loc[:, "ProbeName"].to_list())

        expression_data_collection = []
        expression_indices = []
        sign_expr_data_collection = []
        sign_expr_indices = []
        with gzip.open(self.expr_file, 'rb') as f:
            for i, line in enumerate(f):
                if (i == 0) or (i % self.print_interval == 0):
                    print("\t\tprocessed {} lines\tfound {}/{} expression lines and {}/{} signature genes.".format(i, len(expression_indices), self.eqtl_df.shape[0], len(sign_expr_indices), len(signature_genes)))

                splitted_line = np.array(line.decode().strip('\n').split('\t'))
                index = splitted_line[0]
                data = splitted_line[1:]
                if i == 0:
                    columns = data
                else:
                    hgnc_symbol = None
                    if index in gene_trans_dict.values():
                        hgnc_symbol = gene_trans_dict[index]

                    if hgnc_symbol in expression_interest:
                        expression_indices.append(index)
                        expression_data_collection.append([float(x) for x in data])
                    if hgnc_symbol in signature_genes:
                        sign_expr_indices.append(index)
                        sign_expr_data_collection.append([float(x) for x in data])

        f.close()

        expression_df = pd.DataFrame(expression_data_collection,
                                     index=expression_indices,
                                     columns=columns)
        sign_expr_df = pd.DataFrame(sign_expr_data_collection,
                                    index=sign_expr_indices,
                                    columns=columns)

        return expression_df, sign_expr_df

    def fill_and_reorder(self, df, col):
        left_df = self.eqtl_df.loc[:, col].copy()
        left_df.to_frame(inplace=True)

        combined_df = left_df.merge(df, left_on=col, right_index=True)

        combined_df = combined_df.rename(columns=self.sample_dict)
        combined_df = combined_df[self.sample_order]

        return combined_df

    def clear_variables(self):
        self.geno_file = None
        self.expr_file = None
        self.gene_info_file = None
        self.ens_id = None
        self.hgnc_id = None
        self.force = None

    def get_sign_file(self):
        return self.sign_file

    def get_sign_df(self):
        return self.sign_df

    def get_sign_expr_outpath(self):
        return self.sign_expr_outpath

    def get_gene_info_df(self):
        return self.gene_info_df

    def get_geno_df(self):
        return self.geno_df

    def get_alleles_df(self):
        return self.alleles_df

    def get_expr_df(self):
        return self.expr_df

    def get_sign_expr_df(self):
        return self.sign_expr_df

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > Signature input file: {}".format(self.sign_file))
        print("  > Gene info input file: {}".format(self.gene_info_file))
        print("  > Gene info - ENS column: {}".format(self.ens_id))
        print("  > Gene info - HGNC column: {}".format(self.hgnc_id))
        print("  > eQTL input shape: {}".format(self.eqtl_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Force: {}".format(self.force))
        print("")
