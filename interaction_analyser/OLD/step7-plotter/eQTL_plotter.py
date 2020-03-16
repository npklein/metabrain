#!/usr/bin/env python3

"""
File:         eQTL_plotter.py
Created:      2020/02/26
Last Changed: 2020/03/10
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
from matplotlib.lines import Line2D
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
        self.outdir = os.path.join(os.getcwd(), 'plots_w_corr')
        if group_name is not None:
            self.outdir = os.path.join(self.outdir, group_name)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def start(self, nrows=1):
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

        cov_path = os.path.join(self.indir, 'covariates_table.txt')
        if not os.path.exists(cov_path):
            cov_path = cov_path + '.gz'

        inter_path = os.path.join(self.inter_indir,
                                  'InteractionZScoresMatrix-22Covariates.txt')
        inter_path = "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-12-03-interaction-analyser/analyse_interactions/group_11/output/InteractionZScoresMatrix-40Covariates.txt"

        # if not os.path.exists(inter_path):
        #     inter_path = inter_path + '.gz'
        #
        # # Load the eQTL file.
        # print("Loading eQTL matrix.")
        # eqtl_df = pd.read_csv(eqtl_path, sep="\t", header=0, nrows=nrows)
        # eqtl_df.index = eqtl_df["SNPName"]
        # eqtl_df.index.name = "SNPName"
        # nrows = eqtl_df.shape[0]
        # print(eqtl_df)
        # print("\tShape: {}".format(eqtl_df.shape))
        #
        # # Load the genotype matrix file.
        # print("Loading genotype matrix.")
        # geno_df = pd.read_csv(geno_path, sep="\t", header=0, index_col=0,
        #                       nrows=nrows)
        # geno_df.index.name = "SNPName"
        # print(geno_df)
        # print("\tShape: {}".format(geno_df.shape))
        #
        # # Load the alleles
        # print("Loading alleles matrix.")
        # allele_df = pd.read_csv(allele_path, sep="\t", header=0,
        #                         index_col=0, nrows=nrows)
        # allele_df.index.name = "SNPName"
        # print(allele_df)
        # print("\tShape: {}".format(allele_df.shape))
        #
        # # Load the expression matrix file.
        # print("Loading expression matrix.")
        # expr_df = pd.read_csv(expr_path, sep="\t", header=0, index_col=0,
        #                       nrows=nrows)
        # expr_df.index.name = "SNPName"
        # print(expr_df)
        # print("\tShape: {}".format(expr_df.shape))
        #
        # # Load the covariance matrix file.
        # print("Loading covariance matrix.")
        # cov_df = pd.read_csv(cov_path, sep="\t", header=0, index_col=0)
        # print(cov_df)
        # print("\tShape: {}".format(expr_df.shape))

        # Load the interaction matrix file.
        print("Loading interaction matrix.")
        inter_df = pd.read_csv(inter_path, sep="\t", header=0, index_col=0)
        print(inter_df)
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

        # Calculate the z-score cutoff.
        z_score_cutoff = st.norm.ppf(
            0.05 / (inter_df.shape[0] * inter_df.shape[1]) / 2)

        # Start plotting.

        # Plot heatmap of the interaction matrix.
        print("Creating z-score barplot.")
        self.plot_sum_zscores(inter_df, inter_zscore_outdir)
        print("Creating z-scores distribution plot per covariance.")
        self.plot_distributions(inter_df, z_score_cutoff, inter_zscore_outdir)
        print("Creating z-scores clustermap.")
        self.plot_heatmap(inter_df, inter_zscore_outdir)

        exit()

        # Plot eQTLS.
        group_color_map, value_color_map = self.create_color_map()
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
            allele_map = {0.0: "{}/{}".format(first_allele, first_allele),
                          1.0: "{}/{}".format(first_allele, second_allele),
                          2.0: "{}/{}".format(second_allele, second_allele)}
            data["alleles"] = data["group"].map(allele_map)

            # Add the color.
            data["round_geno"] = data["genotype"].round(2)
            data["value_hue"] = data["round_geno"].map(value_color_map)
            data["group_hue"] = data["group"].map(group_color_map)
            data.drop(["round_geno"], axis=1, inplace=True)

            # Plot a simple eQTL effect.
            self.plot_eqtl_effect(snp_name, probe_name, hgnc_name, eqtl_type,
                                  data, minor_allele, minor_allele_frequency,
                                  first_allele, second_allele,
                                  simple_eqtl_outdir)

            # Check if the SNP has an interaction effect.
            interaction_effect = inter_df.iloc[:, i].to_frame()
            interaction_effect.columns = ["zscore"]
            interaction_effect = interaction_effect.loc[
                                 interaction_effect["zscore"].abs() >= abs(
                                     z_score_cutoff), :]
            interaction_effect = interaction_effect.reindex(
                interaction_effect["zscore"].abs().sort_values(
                    ascending=False).index)

            # Prepare output directory.
            if len(interaction_effect.index) > 0:
                eqtl_interaction_outdir = os.path.join(inter_eqtl_outdir,
                                                       snp_name)
                if not os.path.exists(eqtl_interaction_outdir):
                    os.makedirs(eqtl_interaction_outdir)

                for index, (row,) in interaction_effect.iterrows():
                    eqtl_data = data.copy()
                    cov_data = cov_df.loc[index].to_frame()
                    eqtl_data = eqtl_data.merge(cov_data, left_index=True,
                                                right_index=True)
                    self.plot_interaction_eqtl_effect(snp_name, probe_name,
                                                      hgnc_name,
                                                      eqtl_type, eqtl_data,
                                                      index, row,
                                                      eqtl_interaction_outdir)

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
    def plot_distributions(df, z_score_cutoff, outdir):
        sns.set(style="ticks", color_codes=True)
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
        plt.close()

    @staticmethod
    def plot_sum_zscores(df, outdir, top=10):
        df = df.pow(2)
        sums = df.sum(axis=1).to_frame().reset_index()
        sums.columns = ["index", "counts"]
        sums.sort_values(by=['counts'], ascending=False, inplace=True)

        sns.set()
        fig, ax = plt.subplots(figsize=(11.7, 8.27))
        g = sns.barplot(x="index", y="counts", data=sums, palette="Blues_d")
        g.text(0.5, 1.05,
               'Top Covariates',
               fontsize=16, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.text(0.5, 1.02,
               '',
               fontsize=12, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('sum(z-score^2)',
                     fontsize=12,
                     fontweight='bold')
        g.set_xlabel('covariate',
                     fontsize=12,
                     fontweight='bold')
        ax.tick_params(labelsize=5)
        ax.set_xticks(range(len(sums.index)))
        ax.set_xticklabels(sums["index"], rotation=90)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir, "cov_zscores_barplot.png"))
        plt.close()

        subset = sums.iloc[:top, :]
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
        plt.close()

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
        plt.close()

    @staticmethod
    def plot_eqtl_effect(snp_name, probe_name, hgnc_name, eqtl_type, df,
                         minor_allele, minor_allele_frequency, first_allele,
                         second_allele, outdir):
        """
        """
        # Calculate the correlation.
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            df["genotype"],
            df["expression"])

        # Prepare the figure.
        sns.set()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        fig, ax = plt.subplots()

        # Plot the scatter / box plot.
        sns.regplot(x="genotype", y="expression", data=df,
                    scatter_kws={'facecolors': df['value_hue']},
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
        ax.text(0.5, 1.05,
                '{} {}-eQTL (Pearson r = {:.2f}, '
                'p = {:.2e})'.format(hgnc_name,
                                     eqtl_type,
                                     r_value,
                                     p_value),
                fontsize=16, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        ax.text(0.5, 1.02,
                '',
                fontsize=12, alpha=0.75, ha='center', va='bottom',
                transform=ax.transAxes)
        ax.set_ylabel('{} expression'.format(hgnc_name),
                      fontsize=12,
                      fontweight='bold')
        ax.set_xlabel('SNP {} [minor: {}: {:.2f}'.format(snp_name,
                                                         minor_allele,
                                                         minor_allele_frequency),
                      fontsize=12,
                      fontweight='bold')

        # Show.
        # plt.show()

        # Safe the plot.
        fig.savefig(os.path.join(outdir, "{}_{}_{}.png".format(snp_name,
                                                               probe_name,
                                                               hgnc_name)))
        plt.close()

    @staticmethod
    def plot_interaction_eqtl_effect(snp_name, probe_name, hgnc_name,
                                     eqtl_type, df, cov_name, zscore,
                                     outdir):
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
        fig.savefig(os.path.join(outdir, "{}_{}_{}_{}.png".format(snp_name,
                                                                  probe_name,
                                                                  hgnc_name,
                                                                  cov_name)))
        plt.close()


if __name__ == "__main__":
    # Define main variables.
    # INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
    #                      "output", "2019-11-06-FreezeTwoDotOne",
    #                      "2020-03-03-interaction-analyser",
    #                      "step5-prepare-ia-inputs", "output_p_snp",
    #                      "groups")
    #
    # INTER_INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen",
    #                            "tmp03",
    #                            "output", "2019-11-06-FreezeTwoDotOne",
    #                            "2020-03-03-interaction-analyser",
    #                            "step6-interaction-analyser", "output_wo_corr")

    INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-12-03-interaction-analyser",
                         "matrix_preparation", "output", "create_groups")

    INTER_INDIR = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                               "output", "2019-11-06-FreezeTwoDotOne",
                               "2020-12-03-interaction-analyser",
                               "analyse_interactions", "group_11", "output")

    for GROUP_NAME in next(os.walk(INDIR))[1]:
        if GROUP_NAME != "group_11":
            continue
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
