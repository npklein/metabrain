#!/usr/bin/env python3

"""
File:         eQTL_plotter.py
Created:      2020/02/26
Last Changed: 2020/03/03
Author:       M.Vochteloo

Copyright (C) 2019 M.Vochteloo

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
import gzip

# Local application imports.


# Metadata.
__program__ = "eQTL Plotter"
__author__ = "M. Vochteloo"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class Main:
    """
    Main class of the program.
    """

    def __init__(self, eqtl_file, geno_file, allele_file, expr_file, cov_file,
                 inter_file, interest):
        """
        Initializer method for the main class.

        :param eqtl_file: string, the file containing the eQTL effects of
                          interest
        :param geno_file: string, the genotype data input file.
        :param allele_file: string, the genotype alleles input file.
        :param expr_file: string, the expression data input file.
        :param cov_file: string, the covariate input file.
        :param inter_file: string, the interaction data input file.
        :param interest: string, which SNP's to plot.
        """
        self.eqtl_file = eqtl_file
        self.geno_file = geno_file
        self.allele_file = allele_file
        self.expr_file = expr_file
        self.cov_file = cov_file
        self.inter_file = inter_file
        self.interest = interest
        self.outdir = os.path.join(os.getcwd(), 'output', 'plots')

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        sns.set()

    def start(self, nrows=10):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Load the eQTL file.
        print("Loading eQTL matrix.")
        eqtl_df = pd.read_csv(self.eqtl_file, sep="\t", header=0, nrows=nrows)
        eqtl_df.index = eqtl_df["SNPName"]
        eqtl_df.index.name = "SNPName"
        nrows = eqtl_df.shape[0]
        print("\tShape: {}".format(eqtl_df.shape))

        # Load the genotype matrix file.
        print("Loading genotype matrix.")
        geno_df = pd.read_csv(self.geno_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        geno_df.index.name = "SNPName"
        print("\tShape: {}".format(geno_df.shape))

        # Load the expression matrix file.
        print("Loading expression matrix.")
        expr_df = pd.read_csv(self.expr_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        expr_df.index.name = "SNPName"
        print("\tShape: {}".format(expr_df.shape))

        # Load the interaction matrix file.
        print("Loading interaction matrix.")
        inter_df = pd.read_csv(self.inter_file, sep="\t", header=0, index_col=0,
                               nrows=nrows)
        inter_df.index = [x[:-2] for x in inter_df.index]
        inter_df.index.name = "SNPName"
        print("\tShape: {}".format(expr_df.shape))

        # Load the alleles
        print("Loading alleles matrix.")
        allele_df = pd.read_csv(self.allele_file, sep="\t", header=0,
                                index_col=0, nrows=nrows)
        allele_df.index.name = "SNPName"
        print("\tShape: {}".format(allele_df.shape))

        # Check if shape is identical.
        if (eqtl_df.shape[0] != geno_df.shape[0]) or \
                (eqtl_df.shape[0] != expr_df.shape[0]):
            print("Input matrices rows are not identical length.")
            return
        if geno_df.shape != expr_df.shape:
            print("Genotype and expression matrices are not identical shape.")
            return

        # Check if order is identical.
        if not (eqtl_df.index.identical(geno_df.index)) or \
                not (eqtl_df.index.identical(expr_df.index)) or \
                not (eqtl_df.index.identical(inter_df.index)) or \
                not (eqtl_df.index.identical(allele_df.index)):
            print("Order of SNP's are not identical.")
            return

        # Create the colormap for the simple eqtl plots.
        color_map = self.create_color_map()

        # Prepare output directories.
        simple_eqtl_outdir = os.path.join(self.outdir, "simple_eqtl")
        inter_eqtl_outdir = os.path.join(self.outdir, "interaction_eqtl")
        for dir in [simple_eqtl_outdir, inter_eqtl_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        # Start plotting.
        i = 0
        for index, row in eqtl_df.iterrows():
            if i == "SNPName":
                continue

            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            # Check if gene of interest.
            if (self.interest is not None) and (hgnc_name not in self.interest):
                continue

            # Get the genotype / expression data.
            genotype = geno_df.iloc[i, :].T.to_frame()
            expression = expr_df.iloc[i, :].T.to_frame()
            data = genotype.merge(expression, left_index=True, right_index=True)
            data.columns = ["genotype", "expression"]
            data["group"] = data["genotype"].round(0)

            # Remove missing values.
            data = data.loc[(data['genotype'] >= 0.0) &
                            (data['genotype'] <= 2.0), :]

            # Get the allele data.
            (alleles, minor_allele) = allele_df.iloc[i, :]

            # Determine the genotype order.
            first_allele = snp_name.split(":")[-1].split("_")[0]
            second_allele = snp_name.split(":")[-1].split("_")[1]

            # Check if the major / minor allele gentoypes are correct.
            counts = data["group"].value_counts()
            minor_genotype = counts.idxmin()

            # Determine the minor allele frequency.
            minor_counts = (counts[minor_genotype] * 2) + counts[1.0]
            major_counts = (counts[2.0 - minor_genotype] * 2) + counts[1.0]
            minor_allele_frequency = (minor_counts + major_counts) / (data.shape[0] * 2)

            # Flip the alleles.
            if ((minor_allele == first_allele) and not (
                    minor_genotype == 0.0)) \
                    or ((minor_allele == second_allele) and not (
                    minor_genotype == 2.0)):
                # Flip the genotypes in order to get the genotype labels
                # correct.
                data["genotype"] = 2.0 - data["genotype"]
                data["group"] = 2.0 - data["group"]

            # Check if the SNP has an interaction effect.
            interaction_effect = inter_df.iloc[i, :].to_frame()

            # Plot a simple eQTL effect.
            self.plot_eqtl_effect(snp_name, probe_name, hgnc_name, eqtl_type,
                                  data, minor_allele, minor_allele_frequency,
                                  first_genotype, second_genotype, color_map,
                                  simple_eqtl_outdir)

            # Plot an interaction eQTL effect.
            self.plot_interaction_eqtl_effect(snp_name, probe_name, hgnc_name,
                                              eqtl_type, data,
                                              interaction_effect,
                                              inter_eqtl_outdir)

            i += 1

    @staticmethod
    def create_color_map():
        """
        """
        palette = list(Color("#ABDCA2").range_to(Color("#FFFFFF"), 50)) + \
                  list(Color("#FFFFFF").range_to(Color("#89A3D1"), 50)) + \
                  list(Color("#89A3D1").range_to(Color("#FFFFFF"), 50)) + \
                  list(Color("#FFFFFF").range_to(Color("#F08C72"), 51))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(201))]
        color_map = pd.DataFrame({"value": values,
                                  "hue": colors})
        return color_map

    @staticmethod
    def plot_eqtl_effect(snp_name, probe_name, hgnc_name, eqtl_type, df,
                         minor_allele, minor_allele_frequency, first_allele,
                         second_allele, color_map, outdir):
        """
        """
        # Calculate the correlation.
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            df["genotype"],
            df["expression"])

        # Set the color.
        df["round_geno"] = df["genotype"].round(2)
        df = df.merge(color_map, left_on="round_geno", right_on="value")

        # Prepare the figure.
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter_kws={'facecolors': df['hue']},
                    line_kws={"color": "#D7191C"},
                    ax=ax
                    )
        sns.boxplot(x="group", y="expression", data=df,
                    showfliers=False,
                    zorder=-1,
                    boxprops=dict(alpha=.3),
                    ax=ax)
        plt.setp(ax.artists, edgecolor='k', facecolor='w')
        plt.setp(ax.lines, color='k')

        # Set the other aesthetics.
        ax.set_xticks(range(3))
        ax.set_xticklabels(["{}/{}".format(first_allele, first_allele),
                            "{}/{}".format(first_allele, second_allele),
                            "{}/{}".format(second_allele, second_allele)])
        ax.set_title('{} {}-eQTL (Pearson r = {:.2f}, '
                     'p = {:.2e})'.format(hgnc_name,
                                          eqtl_type,
                                          r_value,
                                          p_value))
        ax.set_ylabel('{} expression'.format(hgnc_name),
                      fontsize=8,
                      fontweight='bold')
        ax.set_xlabel('SNP {} [minor: {}: {:.2f}'.format(snp_name,
                                                         minor_allele,
                                                         minor_allele_frequency),
                      fontsize=8,
                      fontweight='bold')

        # Show.
        # plt.show()

        # Safe the plot.
        fig.savefig(os.path.join(outdir, "{}_{}_{}.png".format(snp_name,
                                                               probe_name,
                                                               hgnc_name)))

    @staticmethod
    def plot_interaction_eqtl_effect(snp_name, probe_name, hgnc_name,
                                     eqtl_type, df, i_effect,
                                     inter_eqtl_outdi):
        print(i_effect)
        exit()


if __name__ == "__main__":
    # Define main variables.
    EQTLS = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output",
                         "2019-11-06-FreezeTwoDotOne",
                         "2020-02-18-eqtls",
                         "cortex-cis-EURandAFR-iterative",
                         "Iteration1",
                         "eQTLProbesFDR0.05-ProbeLevel.txt.gz"
                         )

    GENOTYPE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                            "output", "2019-11-06-FreezeTwoDotOne",
                            "2020-03-03-interaction-analyser",
                            "step4-prepare-matrices", "output",
                            "unmasked", "genotype_table.txt.gz")

    ALLELES = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                           "output", "2019-11-06-FreezeTwoDotOne",
                           "2020-03-03-interaction-analyser",
                           "step4-prepare-matrices", "output",
                           "unmasked", "genotype_alleles.txt.gz")

    EXPRESSION = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output", "2019-11-06-FreezeTwoDotOne",
                              "2020-03-03-interaction-analyser",
                              "step4-prepare-matrices", "output",
                              "unmasked", "expression_table.txt.gz")

    COVARIATE = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                             "output", "2019-11-06-FreezeTwoDotOne",
                             "2020-03-03-interaction-analyser",
                             "step4-prepare-matrices", "output", "unmasked",
                             "covariate_table.txt.gz")

    INTERACTION = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                               "output", "2019-11-06-FreezeTwoDotOne",
                               "2020-02-25-deconvolution", "output",
                               "unmasked", "eQTLInteractionAnalyserOutput",
                               "InteractionZScoresMatrix-0Covariates_T.txt.gz")

    # INTEREST = ["CYB561"]
    INTEREST = None

    # Start the program.
    MAIN = Main(eqtl_file=EQTLS,
                geno_file=GENOTYPE,
                allele_file=ALLELES,
                expr_file=EXPRESSION,
                cov_file=COVARIATE,
                inter_file=INTERACTION,
                interest=INTEREST)
    MAIN.start()
