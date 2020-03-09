#!/usr/bin/env python3

"""
File:         eQTL_plotter.py
Created:      2020/02/26
Last Changed: 2020/03/05
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
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from colour import Color
import scipy.stats as st

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

    def __init__(self, normal_indir, inter_indir, group_name):
        """
        Initializer method for the main class.

        :param normal_indir: string, the input directory.
        :param inter_indir: string, the input directory of the interaction matrix.
        :param group_name, string, the name of the group.
        """
        self.indir = normal_indir
        self.inter_indir = inter_indir
        self.outdir = os.path.join(os.getcwd(), 'plots')
        if group_name is not None:
            self.outdir = os.path.join(self.outdir, group_name)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self, nrows=None):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Interaction input directory: {}".format(self.inter_indir))
        print("  > Output directory: {}".format(self.outdir))
        print("")

        # Construct filenames.
        eqtl_path = os.path.join(self.indir, 'eqtl_table.txt')
        if not os.path.exists(eqtl_path):
            eqtl_path = eqtl_path + '.gz'

        geno_path = os.path.join(self.indir, 'genotype_table.txt')
        if not os.path.exists(geno_path):
            geno_path = geno_path + '.gz'

        allele_path = os.path.join(self.indir, 'genotype_alleles.txt')
        if not os.path.exists(allele_path):
            allele_path = allele_path + '.gz'

        expr_path = os.path.join(self.indir, 'expression_table.txt')
        if not os.path.exists(expr_path):
            expr_path = expr_path + '.gz'

        cov_path = os.path.join(self.indir, 'covariate_table.txt')
        if not os.path.exists(cov_path):
            cov_path = cov_path + '.gz'

        inter_path = os.path.join(self.inter_indir,
                                  'InteractionZScoresMatrix-0Covariates.txt')
        if not os.path.exists(inter_path):
            inter_path = inter_path + '.gz'

        # # Load the eQTL file.
        # print("Loading eQTL matrix.")
        # eqtl_df = pd.read_csv(eqtl_path, sep="\t", header=0, nrows=nrows)
        # eqtl_df.index = eqtl_df["SNPName"]
        # eqtl_df.index.name = "SNPName"
        # nrows = eqtl_df.shape[0]
        # print("\tShape: {}".format(eqtl_df.shape))
        #
        # # Load the genotype matrix file.
        # print("Loading genotype matrix.")
        # geno_df = pd.read_csv(geno_path, sep="\t", header=0, index_col=0,
        #                       nrows=nrows)
        # geno_df.index.name = "SNPName"
        # print("\tShape: {}".format(geno_df.shape))
        #
        # # Load the alleles
        # print("Loading alleles matrix.")
        # allele_df = pd.read_csv(allele_path, sep="\t", header=0,
        #                         index_col=0, nrows=nrows)
        # allele_df.index.name = "SNPName"
        # print("\tShape: {}".format(allele_df.shape))
        #
        # # Load the expression matrix file.
        # print("Loading expression matrix.")
        # expr_df = pd.read_csv(expr_path, sep="\t", header=0, index_col=0,
        #                       nrows=nrows)
        # expr_df.index.name = "SNPName"
        # print("\tShape: {}".format(expr_df.shape))
        #
        # # Load the covariance matrix file.
        # print("Loading covariance matrix.")
        # cov_df = pd.read_csv(cov_path, sep="\t", header=0, index_col=0)
        # print("\tShape: {}".format(expr_df.shape))

        # Load the interaction matrix file.
        print("Loading interaction matrix.")
        inter_df = pd.read_csv(inter_path, sep="\t", header=0, index_col=0)
        print("\tShape: {}".format(inter_df.shape))
        #
        # # Check if shape is identical.
        # if (eqtl_df.shape[0] != geno_df.shape[0]) or \
        #         (eqtl_df.shape[0] != expr_df.shape[0]):
        #     print("Input matrices rows are not identical length.")
        #     return
        # if geno_df.shape != expr_df.shape:
        #     print("Genotype and expression matrices are not identical shape.")
        #     return
        #
        # # Check if SNP order is identical.
        # if not (eqtl_df.index.identical(geno_df.index)) or \
        #         not (eqtl_df.index.identical(expr_df.index)) or \
        #         not (eqtl_df.index.identical(allele_df.index)):
        #     print("Order of SNP's are not identical.")
        #     return
        #
        # # Check if sample order is identical.
        # if not (geno_df.columns.identical(expr_df.columns)) or \
        #         not (geno_df.columns.identical(cov_df.columns)):
        #     print("Order of samples are not identical.")
        #     return

        # Prepare output directories.
        simple_eqtl_outdir = os.path.join(self.outdir, "simple_eqtl")
        inter_eqtl_outdir = os.path.join(self.outdir, "interaction_eqtl")
        inter_zscore_outdir = os.path.join(self.outdir, "interaction_zscores")
        for dir in [simple_eqtl_outdir, inter_eqtl_outdir, inter_zscore_outdir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

        # Start plotting.
        # Plot heatmap of the interaction matrix.
        print("Creating z-score barplot.")
        self.plot_sum_zscores(inter_df, inter_zscore_outdir)
        print("Creating z-scores distribution plot per covariance.")
        self.plot_distributions(inter_df, inter_zscore_outdir)
        print("Creating z-scores clustermap.")
        self.plot_heatmap(inter_df, inter_zscore_outdir)
        return
        # Plot eQTLS.
        color_map = self.create_color_map()
        i = 0
        for index, row in eqtl_df.iterrows():
            if i == "SNPName":
                continue

            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

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
            minor_allele_frequency = (minor_counts + major_counts) / (
                    data.shape[0] * 2)

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
                                  first_allele, second_allele, color_map,
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
    def plot_distributions(df, outdir):
        sns.set(style="ticks", color_codes=True)
        z_score_cutoff = st.norm.ppf(0.05 / (df.shape[0] * df.shape[1]) / 2)
        df = df.T
        dfm = df.melt(var_name='columns')
        g = sns.FacetGrid(dfm, col='columns', col_wrap=10, sharex=False,
                          sharey=False)
        g.map(sns.distplot, 'value')
        g.map(plt.axvline, x=z_score_cutoff, ls='--', c='red')
        g.map(plt.axvline, x=-1 * z_score_cutoff, ls='--', c='red')
        g.set_titles('{col_name}')
        plt.tight_layout()
        g.savefig(os.path.join(outdir, "cov_zscore_distributions.png"))

    @staticmethod
    def plot_sum_zscores(df, outdir, top=10):
        df = df.pow(2)
        sums = df.sum(axis=1).to_frame().reset_index()
        sums.columns = ["index", "counts"]
        sums.sort_values(by=['counts'], ascending=False, inplace=True)

        sns.set()
        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        g = sns.barplot(x="index", y="counts", data=sums, palette="Blues_d")
        g.set_title('Top Covariates')
        g.set_ylabel('sum(z-score^2)',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('covariate',
                     fontsize=8,
                     fontweight='bold')
        ax.tick_params(labelsize=5)
        ax.set_xticks(range(len(sums.index)))
        ax.set_xticklabels(sums["index"], rotation=90)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "cov_zscores_barplot.png"))

        subset = sums.iloc[0:10, :]
        sns.set()
        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        g = sns.barplot(x="index", y="counts", data=subset, palette="Blues_d")
        g.set_title('Top Covariates')
        g.set_ylabel('sum(z-score^2)',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('covariate',
                     fontsize=8,
                     fontweight='bold')
        ax.tick_params(labelsize=10)
        ax.set_xticks(range(len(subset.index)))
        ax.set_xticklabels(subset["index"], rotation=90)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "cov_zscores_barplot_top.png"))

    @staticmethod
    def plot_heatmap(df, outdir):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0, cmap="RdBu_r",
                           yticklabels=True, xticklabels=False,
                           figsize=(12, (.2 * (len(df.index)))))
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=5))
        g.fig.suptitle('Interaction Z-Scores Matrix')
        g.savefig(os.path.join(outdir, "zscores_heatmap.png"))

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
    INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-03-03-interaction-analyser",
                         "step5-prepare-ia-inputs", "output_p_snp",
                         "groups")

    INTER_INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                               "output", "2019-11-06-FreezeTwoDotOne",
                               "2020-03-03-interaction-analyser",
                               "step6-interaction-analyser", "output")

    for GROUP_NAME in next(os.walk(INDIR))[1]:
        GROUP_INDIR = os.path.join(INDIR, GROUP_NAME)
        GROUP_INTER_INDIR = os.path.join(INTER_INDIR, GROUP_NAME)

        # Start the program.
        MAIN = Main(normal_indir=GROUP_INDIR,
                    inter_indir=GROUP_INTER_INDIR,
                    group_name=GROUP_NAME)
        MAIN.start()

    # # Define main variables.
    # INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
    #                      "output", "2019-11-06-FreezeTwoDotOne",
    #                      "2020-03-03-interaction-analyser",
    #                      "step4-prepare-matrices", "output",
    #                      "unmasked")
    #
    # INTER_INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
    #                            "output", "2019-11-06-FreezeTwoDotOne",
    #                            "2020-02-25-deconvolution", "output",
    #                            "unmasked", "eQTLInteractionAnalyserOutput")
    #
    # # Start the program.
    # MAIN = Main(normal_indir=INDIR,
    #             inter_indir=INTER_INDIR,
    #             group_name=None)
    # MAIN.start()
