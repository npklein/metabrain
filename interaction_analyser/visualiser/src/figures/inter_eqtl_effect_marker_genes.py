"""
File:         inter_eqtl_effect_marker_genes.py
Created:      2020/03/17
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
import os

# Third party imports.
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from colour import Color
import scipy.stats as st

# Local application imports.
from src.utilities import prepare_output_dir


class IntereQTLEffectMarkerGenes:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_marker_genes')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        self.eqtl_df = dataset.get_eqtl_df()
        self.geno_df = dataset.get_geno_df()
        self.expr_df = dataset.get_expr_df()
        self.alleles_df = dataset.get_alleles_df()
        self.cov_df = dataset.get_cov_df()
        self.celltypes = dataset.get_celltypes()

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map()

    def start(self):
        print("Plotting interaction eQTL plots for marker genes.")
        self.print_arguments()

        # Determine the columns of the celltypes in the covaraite table.
        marker_indices = []
        for index in self.cov_df.index:
            if ("_" in index) and (
                    index.split("_")[0] in self.celltypes):
                marker_indices.append(index)

        i = 0
        for index, row in self.eqtl_df.iterrows():
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            # Get the genotype / expression data.
            genotype = self.geno_df.iloc[i, :].T.to_frame()
            expression = self.expr_df.iloc[i, :].T.to_frame()
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]
            data["group"] = data["genotype"].round(0)

            # Remove missing values.
            data = data.loc[(data['genotype'] >= 0.0) &
                            (data['genotype'] <= 2.0), :]

            # Get the allele data.
            (alleles, minor_allele) = self.alleles_df.iloc[i, :]

            # Determine the genotype order.
            first_allele = snp_name.split(":")[-1].split("_")[0]
            second_allele = snp_name.split(":")[-1].split("_")[1]

            # Check if the major / minor allele gentoypes are correct.
            counts = data["group"].value_counts()
            minor_genotype = counts.idxmin()

            # Flip the alleles.
            if ((minor_allele == first_allele) and not (
                    minor_genotype == 0.0)) \
                    or ((minor_allele == second_allele) and not (
                    minor_genotype == 2.0)):
                # Flip the genotypes in order to get the genotype labels
                # correct.
                data["genotype"] = 2.0 - data["genotype"]
                data["group"] = 2.0 - data["group"]
            allele_map = {0.0: "{}/{}".format(first_allele, first_allele),
                          1.0: "{}/{}".format(first_allele, second_allele),
                          2.0: "{}/{}".format(second_allele, second_allele)}
            data["alleles"] = data["group"].map(allele_map)

            # Add the color.
            data["round_geno"] = data["genotype"].round(2)
            data["value_hue"] = data["round_geno"].map(self.value_color_map)
            data["group_hue"] = data["group"].map(self.group_color_map)
            data.drop(["round_geno"], axis=1, inplace=True)

            # Get the covariates of the marker genes.
            covariates = self.cov_df.loc[marker_indices, :]

            self.plot(snp_name, probe_name, hgnc_name, eqtl_type,
                      data, covariates, self.celltypes, self.outdir)

            i += 1

    @staticmethod
    def create_color_map():
        """
        """
        palette = list(Color("#ABDCA2").range_to(Color("#CBE9C5"), 50)) + \
                  list(Color("#B1C2E1").range_to(Color("#89A3D1"), 50)) + \
                  list(Color("#89A3D1").range_to(Color("#B1C2E1"), 50)) + \
                  list(Color("#F6BCAD").range_to(Color("#F08C72"), 51))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(201))]
        group_color_map = {0.0: "#ABDCA2", 1.0: "#89A3D1", 2.0: "#F08C72"}
        value_color_map = {}
        for val, col in zip(values, colors):
            value_color_map[val] = col
        return group_color_map, value_color_map

    @staticmethod
    def plot(snp_name, probe_name, hgnc_name, eqtl_type, df, markers,
             celltypes, outdir):
        """
        """
        print(celltypes)
        print(df)
        print(markers)

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_df.shape))
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
