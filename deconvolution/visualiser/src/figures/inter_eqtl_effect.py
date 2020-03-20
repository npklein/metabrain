"""
File:         inter_eqtl_effect.py
Created:      2020/03/16
Last Changed: 2020/03/20
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
from colour import Color
import os

# Third party imports.
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLEffect:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.eqtl_df = dataset.get_eqtl_df()
        self.geno_df = dataset.get_geno_df()
        self.expr_df = dataset.get_expr_df()
        self.alleles_df = dataset.get_alleles_df()
        self.cov_df = dataset.get_cov_df()
        self.inter_df = dataset.get_inter_df()

        # Create color map.
        self.group_color_map, self.value_color_map = self.create_color_map()

    def start(self):
        print("Plotting interaction eQTL plots.")
        self.print_arguments()

        # Calculate the z-score cutoff.
        z_score_cutoff = stats.norm.ppf(
            0.05 / (self.inter_df.shape[0] * self.inter_df.shape[1]) / 2)

        for i, (index, row) in enumerate(self.eqtl_df.iterrows()):
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            print("\tWorking on: {}\t{}\t{} [{}/{} "
                  "{:.2f}%]".format(snp_name, probe_name, hgnc_name,
                                    i + 1,
                                    self.eqtl_df.shape[0],
                                    (100 / self.eqtl_df.shape[0]) * (i + 1)))

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

            # Check if the SNP has an interaction effect.
            interaction_effect = self.inter_df.iloc[:, i].to_frame()
            interaction_effect.columns = ["zscore"]
            interaction_effect = interaction_effect.loc[
                                 interaction_effect["zscore"].abs() >= abs(
                                     z_score_cutoff), :]
            interaction_effect = interaction_effect.reindex(
                interaction_effect["zscore"].abs().sort_values(
                    ascending=False).index)

            # Prepare output directory.
            if len(interaction_effect.index) > 0:
                eqtl_interaction_outdir = os.path.join(self.outdir,
                                                       "{}_{}".format(i,
                                                                      snp_name))
                if not os.path.exists(eqtl_interaction_outdir):
                    os.makedirs(eqtl_interaction_outdir)

                count = 0
                for index, (row,) in interaction_effect.iterrows():
                    eqtl_data = data.copy()
                    cov_data = self.cov_df.loc[index].to_frame()
                    eqtl_data = eqtl_data.merge(cov_data, left_index=True,
                                                right_index=True)
                    self.plot(snp_name, probe_name, hgnc_name, eqtl_type,
                              eqtl_data, index, row, count,
                              eqtl_interaction_outdir)
                    count += 1

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
    def plot(snp_name, probe_name, hgnc_name, eqtl_type, df, cov_name, zscore,
             count, outdir):
        """
        """
        # calculate axis limits.
        min = df["expression"].min() * 1.1
        max = df["expression"].max() * 1.5

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sns.despine(fig=fig, ax=ax)

        for i, allele in enumerate(df["alleles"].unique()):
            # Calculate the correlation.
            subset = df.loc[df["alleles"] == allele, :].copy()
            coef, p = stats.spearmanr(subset["expression"], subset[cov_name])
            color = subset["group_hue"][0]

            # Plot.
            sns.regplot(x=cov_name, y="expression", data=subset,
                        scatter_kws={'facecolors': subset['value_hue'],
                                     'edgecolors': subset['group_hue']},
                        line_kws={"color": color},
                        ax=ax
                        )

            # Add the text.
            ax.set(ylim=(min, max))
            ax.annotate(
                '{}: r = {:.2e}, p = {:.2e}'.format(allele, coef, p),
                xy=(0.03, 0.94 - ((i / 100) * 3)),
                xycoords=ax.transAxes,
                color=color,
                fontsize=12,
                fontweight='bold')

        ax.text(0.5, 1.05,
                '{} {}-eQTL Interaction with {} '
                '[z-score: {:.2f}]'.format(hgnc_name,
                                           eqtl_type,
                                           cov_name,
                                           zscore),
                fontsize=16, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                'SNPName: {}    ProbeName:{}'.format(snp_name, probe_name),
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)

        ax.set_ylabel('{} ({}) expression'.format(probe_name, hgnc_name),
                      fontsize=12,
                      fontweight='bold')
        ax.set_xlabel(cov_name,
                      fontsize=12,
                      fontweight='bold')

        # Safe the plot.
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_inter_eqtl_{}_{}_{}_{}.png".format(
                                     count,
                                     snp_name,
                                     probe_name,
                                     hgnc_name,
                                     cov_name)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_df.shape))
        print("  > Genotype matrix shape: {}".format(self.geno_df.shape))
        print("  > Alleles matrix shape: {}".format(self.alleles_df.shape))
        print("  > Expression matrix shape: {}".format(self.expr_df.shape))
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
