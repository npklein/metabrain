#!/usr/bin/env python3

"""
File:         prepare_ia_inputs_cc_per_snp.py
Created:      2020/03/05
Last Changed: 2020/03/10
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
from itertools import groupby
import os

# Third party imports.
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.


# Metadata.
__program__ = "Prepare Interaction Analyser Inputs (Complete Case per SNP)"
__author__ = "M. Vochteloo"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class Group:
    def __init__(self, id, samples, cohorts_list):
        self.id = id
        self.samples = samples
        self.cohort_counts = self.calc_cohort_counts(cohorts_list)

        self.snp_indices = np.array([], dtype=np.int16)
        self.sample_indices = None
        self.eqtls = []

    @staticmethod
    def calc_cohort_counts(cohort_list):
        counts = {}
        for key, group in groupby(cohort_list):
            counts[key] = len(list(group))
        return counts

    def add_eqtl(self, eqtl):
        self.eqtls.append(eqtl)
        self.add_snp_index(eqtl.get_snp_index())
        self.set_sample_indices(eqtl.get_sample_indices())

    def add_snp_index(self, index):
        self.snp_indices = np.append(self.snp_indices, index)

    def set_sample_indices(self, indices):
        if self.sample_indices is None:
            self.sample_indices = indices

    def get_id(self):
        return self.id

    def get_samples(self):
        return self.samples

    def get_cohort_counts(self):
        return self.cohort_counts

    def get_snp_indices(self):
        return self.snp_indices

    def get_sample_indices(self):
        return self.sample_indices

    def get_eqtls(self):
        return self.eqtls

    def get_n_eqtls(self):
        return len(self.eqtls)

    def get_n_samples(self):
        return len(self.sample_indices)

    def matches(self, sample_indices):
        return np.array_equal(self.sample_indices, sample_indices)


class EQTL:
    def __init__(self, name, snp_index, genotype, expression):
        self.name = name
        self.snp_index = snp_index

        df = self.combine_data(genotype, expression)
        self.samples, self.sample_indices = self.filter_data(df)
        self.n_missing = genotype.size - len(self.samples)
        del df

        self.group = None

    @staticmethod
    def combine_data(genotype, expression):
        genotype_df = genotype.T.to_frame()
        expression_df = expression.T.to_frame()
        data = genotype_df.merge(expression_df,
                                 left_index=True,
                                 right_index=True)
        data.columns = ["genotype", "expression"]
        del genotype_df, expression_df
        return data

    @staticmethod
    def filter_data(df):
        tmp_df = df.copy()
        tmp_df.replace(-1, np.nan, inplace=True)
        indices = np.arange(tmp_df.shape[0])
        sample_indices = indices[~tmp_df.isnull().any(axis=1).values]
        tmp_df.dropna(inplace=True)
        samples = tmp_df.index.to_list()
        del tmp_df
        return samples, sample_indices

    def get_name(self):
        return self.name

    def get_snp_index(self):
        return self.snp_index

    def get_samples(self):
        return self.samples

    def get_sample_indices(self):
        return self.sample_indices

    def get_n_samples(self):
        return len(self.samples)

    def get_n_missing(self):
        return self.n_missing


class Main:
    """
    Main class of the program.
    """

    def __init__(self, geno_file, allele_file, expr_file, cov_file, cohorts,
                 eqtl_file, eqtl_ia):
        """
        Initializer method for the main class.

        :param geno_file: string, the genotype data input file.
        :param allele_file: string, the alleles data input file.
        :param expr_file: string, the expression data input file.
        :param cov_file: string, the covariates data input file.
        :param cohorts: list, the unique cohorts in the data.
        :param eqtl_file: string, the eQTL data input file.
        :param eqtl_ia: string, the path to the interaction analyser.
        """
        self.geno_file = geno_file
        self.allele_file = allele_file
        self.expr_file = expr_file
        self.cov_file = cov_file
        self.cohorts = cohorts
        self.eqtl_file = eqtl_file
        self.eqtl_ia = eqtl_ia

        self.outdir = os.path.join(os.getcwd(), 'output_p_snp')
        self.plot_dir = os.path.join(self.outdir, 'plots')
        self.group_dir = os.path.join(self.outdir, 'groups')

        for dir in [self.outdir, self.plot_dir, self.group_dir]:
            if not os.path.exists(dir):
                os.makedirs(dir)

    def start(self, nrows=None, min_eqtls_in_group=1000):
        """
        Main method for the main class. Does all the work.

        :param nrows: int, the number of rows to parse of the input file.
                      used for development.
        :param min_eqtls_in_group: int, the minimal number of eQTLs in the
                                    group.
        """
        # Print arguments.
        print("Arguments:")
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > Covariate input file: {}".format(self.cov_file))
        print("  > Cohorts: {}".format(self.cohorts))
        print("  > EQTL input file: {}".format(self.eqtl_file))
        print("  > eQTL Interaction Analyser: {}".format(self.eqtl_ia))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Min. number of eQTLs in a group: {}".format(
            min_eqtls_in_group))
        print("")

        # Load the genotype data.
        print("Loading genotype matrix.")
        geno_df = pd.read_csv(self.geno_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        print("\tShape: {}".format(geno_df.shape))

        # Load the expression data.
        print("Loading expression matrix.")
        expr_df = pd.read_csv(self.expr_file, sep="\t", header=0, index_col=0,
                              nrows=nrows)
        print("\tShape: {}".format(expr_df.shape))

        # Load the coviarate data.
        print("Loading covariate matrix.")
        cov_df = pd.read_csv(self.cov_file, sep="\t", header=0, index_col=0,
                             nrows=nrows)
        print("\tShape: {}".format(cov_df.shape))

        # Create a sample-cohort dict: {sample: cohort, etc.}
        translate_dict = {}
        cohort_counts = {}
        for cohort in self.cohorts:
            samples_in_cohort = cov_df.loc[
                cohort, (cov_df.loc[cohort] == 1) | (cov_df.loc[cohort] == "1")].index.to_list()
            cohort_counts[cohort] = len(samples_in_cohort)
            for sample in samples_in_cohort:
                translate_dict[sample] = cohort

        # Check if the shape is identical.
        if geno_df.shape != expr_df.shape:
            print("Genotype and expression matrices are not identical shape.")
            return

        # Loop over the eQTL.
        print("Processing.")
        n_eqtls = geno_df.shape[0]
        groups = []
        new_group_id = 0
        for i in range(n_eqtls):
            if (i % 250 == 0) or (i == (n_eqtls - 1)):
                print("\t Processing {}/{} [{:.2f}%]".format(i, n_eqtls - 1,
                                                             (100 / (
                                                                     n_eqtls - 1)) * i))

            # Get the snp name.
            snp_name = geno_df.index[i]

            # Get the genotype / expression data.
            genotype = geno_df.iloc[i, :].copy()
            expression = expr_df.iloc[i, :].copy()

            # Create an eQTL object.
            new_eqtl = EQTL(snp_name, i, genotype, expression)

            # Get the samples indices of the eQTl.
            samples = new_eqtl.get_samples()
            samples_indices = new_eqtl.get_sample_indices()

            # Assign the group.
            matches = False
            if groups:
                # Check if there is a group with these samples.
                for group in groups:
                    if group.matches(samples_indices):
                        group.add_eqtl(new_eqtl)
                        matches = True
                        break

            # Add a new group.
            if not matches:
                cohorts_list = []
                for sample in samples:
                    if sample in translate_dict.keys():
                        cohorts_list.append(translate_dict[sample])
                new_group = Group(new_group_id, samples, cohorts_list)
                new_group.add_eqtl(new_eqtl)
                groups.append(new_group)
                new_group_id = new_group_id + 1

        # Plot counts per sample.
        sample_counts = pd.DataFrame(0,
                                     index=geno_df.index,
                                     columns=['present', 'missing', 'complete'])
        sample_counts.index.name = "SNPName"
        n_groups = 0
        n_eqtls = 0
        for group in groups:
            n_groups += 1
            for eqtl in group.get_eqtls():
                n_eqtls += 1
                sample_counts.at[
                    eqtl.get_name(), 'present'] = eqtl.get_n_samples()
                sample_counts.at[
                    eqtl.get_name(), 'missing'] = eqtl.get_n_missing()
                sample_counts.at[eqtl.get_name(), 'complete'] = geno_df.shape[1]
        print("Number of groups: {}".format(n_groups))
        print("Number of eQTLs: {}".format(n_eqtls))
        sample_counts.to_csv(os.path.join(self.outdir, 'sample_counts.txt.gz'),
                             sep="\t",
                             compression="gzip")
        self.plot_distribution(sample_counts)

        # Create group dataframes.
        group_counts = pd.DataFrame(0,
                                    index=['total'] + list(range(new_group_id)),
                                    columns=['n_eqtls', 'n_samples'] + self.cohorts)
        group_counts.index.name = "GroupName"
        for cohort, count in cohort_counts.items():
            group_counts.at['total', cohort] = count
        group_counts.at['total', 'n_eqtls'] = geno_df.shape[0]
        group_counts.at['total', 'n_samples'] = geno_df.shape[1]
        for group in groups:
            group_counts.at[group.get_id(), 'n_eqtls'] = group.get_n_eqtls()
            group_counts.at[group.get_id(), 'n_samples'] = group.get_n_samples()
            for cohort, count in group.get_cohort_counts().items():
                group_counts.at[group.get_id(), cohort] = count
        # group_counts = group_counts.loc[group_counts['n_eqtls'] > 1, :]
        group_counts.sort_values(by='n_eqtls', ascending=False, inplace=True)
        group_counts.to_csv(os.path.join(self.group_dir, 'group_counts.txt.gz'),
                            sep="\t",
                            compression="gzip")
        print("Groups:")
        print(group_counts.iloc[0:15, :])
        subset = group_counts.loc[group_counts['n_eqtls'] >= min_eqtls_in_group,
                 :]
        self.plot_stacked_bar_graph(subset)

        # Load the eqtl data.
        print("\tLoading covariate matrix.")
        eqtl_df = pd.read_csv(self.eqtl_file, sep="\t", header=0,
                              nrows=nrows)
        print("\t\tShape: {}".format(eqtl_df.shape))

        # Load the alleles data.
        print("\tLoading alleles matrix.")
        allele_df = pd.read_csv(self.allele_file, sep="\t", header=0,
                                nrows=nrows)
        print("\t\tShape: {}".format(eqtl_df.shape))

        # Safe output files.
        for group in groups:
            if group.get_n_eqtls() >= min_eqtls_in_group:
                print("Saving files for group: {}.".format(group.get_id()))
                print("Group has {} eQTL's for {} samples.".format(
                    group.get_n_eqtls(),
                    group.get_n_samples()))
                # Define and prepare the output directory.
                group_dir = os.path.join(self.group_dir,
                                         'group_{}'.format(group.get_id()))
                if not os.path.exists(group_dir):
                    os.makedirs(group_dir)

                # Get the group indices.
                snp_mask = group.get_snp_indices()
                sample_mask = group.get_sample_indices()

                # Write the uncompressed genotype file.
                group_geno = geno_df.iloc[snp_mask, sample_mask].copy()
                print("\tGenotype matrix shape:   {}".format(group_geno.shape))
                # validate that there are no -1 left.
                tmp_df = group_geno.iloc[:, 1:].copy()
                tmp_df.replace(-1, np.nan, inplace=True)
                if np.count_nonzero(np.isnan(tmp_df)) != 0:
                    print("Matrix still has missing genotypes.")
                    exit()
                del tmp_df
                group_geno_path = os.path.join(group_dir, 'genotype_table.txt')
                group_geno.to_csv(group_geno_path, sep="\t")
                del group_geno

                # Write the uncompressed expression file.
                group_expr = expr_df.iloc[snp_mask, sample_mask].copy()
                print("\tExpression matrix shape: {}".format(group_expr.shape))
                group_expr_path = os.path.join(group_dir,
                                               'expression_table.txt')
                group_expr.to_csv(group_expr_path, sep="\t")
                del group_expr

                # Safe the covariance table of the group.
                group_cov = cov_df.iloc[:, sample_mask].copy()
                group_cov = group_cov.loc[group_cov.std(axis=1) != 0, :]
                print("\tCovariate matrix shape:  {}".format(group_cov.shape))
                group_cov_path = os.path.join(group_dir, 'covariate_table.txt')
                group_cov.to_csv(group_cov_path, sep="\t")
                del group_cov

                # Safe the eqtl table of the group.
                group_eqtl = eqtl_df.iloc[snp_mask, :].copy()
                print("\tEQTL matrix shape:       {}".format(group_eqtl.shape))
                group_eqtl.to_csv(os.path.join(group_dir, 'eqtl_table.txt.gz'),
                                  index=False, sep="\t", compression="gzip")
                del group_eqtl

                # Safe the allele table of the group.
                group_alleles = allele_df.iloc[snp_mask, :].copy()
                print("\tAlleles matrix shape:       {}".format(group_alleles.shape))
                group_alleles.to_csv(os.path.join(group_dir, 'genotype_alleles.txt.gz'),
                                     index=False, sep="\t", compression="gzip")
                del group_alleles

                # Convert genotype / expression / covariance files to binary.
                self.create_binary_file(group_geno_path, 'Genotypes')
                self.create_binary_file(group_expr_path, 'Expression')
                self.create_binary_file(group_cov_path, 'Covariates')

    def plot_distribution(self, df):
        print("Plotting distribution.")
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        g = sns.distplot(df[['present']], ax=ax, bins=25)
        g.set_title('Number of eQTLs per sample size')
        g.set_ylabel('Frequency',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('Sample size',
                     fontsize=8,
                     fontweight='bold')
        fig.savefig(os.path.join(self.plot_dir,
                                 "sample_distribution.png"))

    def plot_stacked_bar_graph(self, df):
        n_groups = df.shape[0]
        group_ids = df.index.to_list()
        eqtl_counts = df.pop('n_eqtls').to_list()
        sample_couns = df.pop('n_samples').to_list()
        molten_df = pd.melt(df.reset_index(), id_vars='GroupName')

        print("Plotting stacked bar graph.")
        fig, ax = plt.subplots()
        sns.set_style("darkgrid", {"axes.facecolor": ".9"})
        g = sns.barplot(x="GroupName", y="value", hue="variable",
                        data=molten_df, palette="deep")
        g.set_title('Grouping of eQTLs')
        g.set_ylabel('Sample size',
                     fontsize=8,
                     fontweight='bold')
        g.set_xlabel('Group',
                     fontsize=8,
                     fontweight='bold')
        ax.set_xticks(range(n_groups))
        ax.set_xticklabels(["Group: {}\nN-eQTLs: {}\n"
                            "N-Samples: {}".format(group_ids[i],
                                                   eqtl_counts[i],
                                                   sample_couns[i])
                            for i in range(n_groups)])
        lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        fig.savefig(os.path.join(self.plot_dir,
                                 "group_info.png"),
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

    def create_binary_file(self, inpath, out_filename):
        """
        Method for creating a binary file for the Interaction Analyser input.

        :param inpath: string, the input file path.
        :param out_filename: string, the output filename.
        :return:
        """
        outpath = os.path.join(os.path.dirname(inpath),
                               out_filename + '.binary')

        # Convert to binary.
        command = 'java -jar {} --convertMatrix -i {} -o {}'.format(
            self.eqtl_ia, inpath, outpath)
        print("\t{}".format(command))
        os.system(command)


if __name__ == "__main__":
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

    COVARIATES = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                              "output", "2019-11-06-FreezeTwoDotOne",
                              "2020-03-03-interaction-analyser",
                              "step4-prepare-matrices", "output",
                              "unmasked", "covariate_table.txt.gz")

    COHORTS = ["AMPAD-MSBB-V2-AFR",
               "CMC-AFR",
               "LIBD_1M-AFR",
               "LIBD_h650-AFR",
               "AMPAD-MAYO-V2-EUR",
               "AMPAD-MSBB-V2-EUR",
               "AMPAD-ROSMAP-V2-EUR",
               "BrainGVEX-V2-EUR",
               "CMC-EUR",
               "GTEx-EUR",
               "GVEx",
               "LIBD_1M-EUR",
               "LIBD_h650-EUR",
               "NABEC-H550-EUR",
               "NABEC-H610-EUR",
               "TargetALS-EUR",
               "UCLA_ASD-EUR",
               "ENA-EU"]

    EQTLS = os.path.join(os.path.sep, "groups", "umcg-biogen", "tmp03",
                         "output", "2019-11-06-FreezeTwoDotOne",
                         "2020-03-03-interaction-analyser",
                         "step1-combine-eQTLprobe-files", "output",
                         "eQTLProbesFDR0.05-ProbeLevel_combined.txt.gz")

    EQTL_IA = os.path.join(os.path.sep, "groups", "umcg-biogen",
                           "tmp03", "output", "2019-11-06-FreezeTwoDotOne",
                           "2020-03-03-interaction-analyser",
                           "step6-interaction-analyser",
                           "eQTLInteractionAnalyser-1.2-SNAPSHOT-jar-with-dependencies.jar")

    # Start the program.
    MAIN = Main(geno_file=GENOTYPE,
                allele_file=ALLELES,
                expr_file=EXPRESSION,
                cov_file=COVARIATES,
                cohorts=COHORTS,
                eqtl_file=EQTLS,
                eqtl_ia=EQTL_IA)

    MAIN.start()
