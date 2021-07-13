#!/usr/bin/env python3

"""
File:         create_decon_eqtl_permutation_jobs.py
Created:      2021/06/04
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
import random
import os

# Third party imports.
import numpy as np
import pandas as pd
from statsmodels.stats import multitest

# Local application imports.

"""
Syntax:
./create_decon_eqtl_permutation_jobs.py -de /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -eq /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ge /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -ex /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/cis/cortex/expression_table/2020-07-16-MetaBrainDeconQtlGenes.TMM.SampSelect.ZeroVarRemov.covRemoved.expAdded.txt -cc /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt.gz -co /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/matrix_preparation/cortex_eur_cis/create_cohort_matrix/cohort_matrix.txt.gz
"""

# Metadata
__program__ = "Create Decon-eQTL Permutation Jobs"
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


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.decon_path = getattr(arguments, 'decon')
        self.eqtl_path = getattr(arguments, 'eqtl')
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.cc_path = getattr(arguments, 'cell_counts')
        self.cohort_path = getattr(arguments, 'cohort')
        self.alpha = getattr(arguments, 'alpha')
        self.n_permutations = getattr(arguments, 'permutations')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent.absolute()), "decon_eqtl_permutation_per_cohort")
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
        parser.add_argument("-de",
                            "--decon",
                            type=str,
                            required=True,
                            help="The path to the decon output "
                                 "matrix")
        parser.add_argument("-eq",
                            "--eqtl",
                            type=str,
                            required=True,
                            help="The path to the eqtl matrix")
        parser.add_argument("-ge",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-ex",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-cc",
                            "--cell_counts",
                            type=str,
                            required=True,
                            help="The path to the cell counts "
                                 "matrix")
        parser.add_argument("-co",
                            "--cohort",
                            type=str,
                            required=True,
                            help="The path to the cohort matrix.")
        parser.add_argument("-a",
                            "--alpha",
                            type=float,
                            required=False,
                            default=0.05,
                            help="The significance cut-off. Default: 0.05.")
        parser.add_argument("-p",
                            "--permutations",
                            type=int,
                            default=1000,
                            help="The number of permutations to run.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("Prepare output directories.")
        decon_outdir = os.path.join(self.outdir, "decon_output")
        real_decon_outdir = os.path.join(decon_outdir, "real")
        perm_genotype_outdir = os.path.join(self.outdir, "perm_genotypes")
        jobsdir = os.path.join(self.outdir, "jobs")
        jobs_outdir = os.path.join(self.outdir, "jobs", "output")
        for outdir in [decon_outdir, real_decon_outdir, perm_genotype_outdir, jobsdir, jobs_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        print("Loading eQTL and Decon-QTL data.")
        eqtl_df = self.load_file(self.eqtl_path, header=0, index_col=None)
        decon_df = self.load_file(self.decon_path, header=0, index_col=0)

        print("Calculating FDR.")
        decon_fdr_df, decon_columns = self.bh_correct(decon_df)
        del decon_df

        print("Merging eQTL data.")
        eqtl_df = eqtl_df.merge(decon_fdr_df, on=["ProbeName", "SNPName"])
        del decon_fdr_df

        print("Finding indices with significant hit.")
        mask = eqtl_df[decon_columns].min(axis=1) < self.alpha
        eqtl_df = eqtl_df.loc[mask, :]
        n_eqtls = np.sum(mask)
        n_ieqtls = np.sum(eqtl_df[decon_columns] < self.alpha).sum()
        print("\t{} eQTLs have a total of {} interaction.".format(n_eqtls,
                                                                  n_ieqtls))

        print("Loading other data")
        skiprows = None
        if np.sum(mask) > 0:
            skiprows = [x for x in np.arange(1, np.size(mask) + 1)[~mask]]
        geno_df = self.load_file(self.geno_path, header=0, index_col=0, skiprows=skiprows, nrows=np.size(mask) + 1)
        expr_df = self.load_file(self.expr_path, header=0, index_col=0, skiprows=skiprows, nrows=np.size(mask) + 1)
        cc_df = self.load_file(self.cc_path, header=0, index_col=0)
        cohort_df = self.load_file(self.cohort_path, header=0, index_col=0)
        cohort_df = cohort_df.loc[cohort_df.sum(axis=1) > 0, :]

        print("\tCohorts:")
        for index, row in cohort_df.iterrows():
            print("\t{}: n-samples: {}".format(index, cohort_df.loc[index, :].sum()))

        # select samples from expr matrix.
        # this is done manually cause the ZeroVarRemov.covRemoved.expAdded.txt
        # file that decon needs suddenly has a weird added column.
        expr_df = expr_df.loc[:, geno_df.columns.tolist()]

        print("Validating input.")
        self.validate_data(eqtl_df=eqtl_df,
                           geno_df=geno_df,
                           expr_df=expr_df,
                           cc_df=cc_df,
                           cohort_df=cohort_df
                           )

        print("Saving uncompressed data")
        snp_gene_df = eqtl_df.loc[mask, ["ProbeName", "SNPName"]]
        snp_gene_df.columns = ["ProbeName", "SNP"]
        snp_gene_filepath = os.path.join(self.outdir, "snp_gene_list.txt")
        self.save_file(snp_gene_df, outpath=snp_gene_filepath, index=False)

        geno_filepath = os.path.join(self.outdir, "genotype.txt")
        self.save_file(geno_df, outpath=geno_filepath)

        expr_filepath = os.path.join(self.outdir, "expression.txt")
        self.save_file(expr_df, outpath=expr_filepath)

        cc_filepath = os.path.join(self.outdir, "cell_counts.txt")
        self.save_file(cc_df, outpath=cc_filepath)
        del eqtl_df, snp_gene_df, expr_df, cc_df

        print("Saving the primary job file.")
        self.create_job_file(job_name="real",
                             jobsdir=jobsdir,
                             jobs_outdir=jobs_outdir,
                             results_outdir=real_decon_outdir,
                             cellcount=cc_filepath,
                             expression=expr_filepath,
                             snps_to_test=snp_gene_filepath,
                             genotype=geno_filepath,
                             )

        print("Create permutations")
        perm_order_m = self.create_perm_orders(n_permutations=self.n_permutations,
                                               cohort_df=cohort_df)
        print(pd.DataFrame(perm_order_m))
        self.save_file(pd.DataFrame(perm_order_m),
                       outpath=os.path.join(self.outdir, "perm_order.txt.gz"),
                       index=False, header=False)

        print("Saving permutations genotype files.")
        digit_length = len(str(self.n_permutations))
        geno_columns = geno_df.columns.tolist()
        for i in range(perm_order_m.shape[0]):
            perm_name = "perm{:0{}d}".format(i, digit_length)

            perm_decon_outdir = os.path.join(decon_outdir, perm_name)
            if not os.path.exists(perm_decon_outdir):
                os.makedirs(perm_decon_outdir)

            # Shuffle the genotype.
            perm_geno_df = geno_df.iloc[:, perm_order_m[i, :]]
            perm_geno_df.columns = geno_columns
            perm_geno_filepath = os.path.join(perm_genotype_outdir, "genotype_{}.txt".format(perm_name))
            self.save_file(perm_geno_df, outpath=perm_geno_filepath)

            self.create_job_file(job_name=perm_name,
                                 jobsdir=jobsdir,
                                 jobs_outdir=jobs_outdir,
                                 results_outdir=perm_decon_outdir,
                                 cellcount=cc_filepath,
                                 expression=expr_filepath,
                                 snps_to_test=snp_gene_filepath,
                                 genotype=perm_geno_filepath,
                                 )

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
    def bh_correct(pvalue_df):
        df = pvalue_df.copy()
        fdr_data = []
        indices = []
        for col in df.columns:
            if col.endswith("_pvalue"):
                fdr_data.append(
                    multitest.multipletests(df.loc[:, col], method='fdr_bh')[1])
                indices.append(col.replace("_pvalue", ""))
        fdr_df = pd.DataFrame(fdr_data, index=indices, columns=df.index).T
        columns = fdr_df.columns.tolist()

        # Split the index.
        probe_names = []
        snp_names = []
        for index in fdr_df.index:
            probe_names.append(index.split("_")[0])
            snp_names.append("_".join(index.split("_")[1:]))
        fdr_df.insert(0, "ProbeName", probe_names)
        fdr_df.insert(1, "SNPName", snp_names)

        return fdr_df, columns

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

    @staticmethod
    def validate_data(eqtl_df, geno_df, expr_df, cc_df, cohort_df):
        snp_reference = eqtl_df["SNPName"].values.tolist()
        probe_reference = eqtl_df["ProbeName"].values.tolist()

        if not geno_df.index.tolist() == snp_reference:
            print("The genotype file indices do not match the "
                           "eQTL file.")
            exit()

        if not expr_df.index.tolist() == probe_reference:
            print("The expression file indices do not match the "
                  "eQTL file.")
            exit()

        if not geno_df.columns.tolist() == expr_df.columns.tolist():
            print("The genotype file header does not match the "
                  "expression file header.")
            exit()

        if not geno_df.columns.tolist() == cc_df.index.tolist():
            print("The genotype file header does not match the "
                  "cell count file index.")
            exit()

        if not geno_df.columns.tolist() == cohort_df.columns.tolist():
            print("The genotype file header does not match the "
                  "cohort file header.")
            exit()

        if not cc_df.columns.tolist() == eqtl_df.columns[-cc_df.shape[1]:].tolist():
            print("The cell count file header does not match the "
                  "deconvolution file header.")
            exit()

        if not cohort_df.sum(axis=0).sum() == cohort_df.shape[1]:
            print("The cohort file is invalid.")
            exit()

    @staticmethod
    def create_job_file(job_name, jobsdir, jobs_outdir, results_outdir,
                        cellcount, expression, snps_to_test, genotype,
                        ParallelGCThreads=2, initial_heap_size=1000,
                        max_heap_size=1500, base_name="decon_cis_cortex_",
                        time="00:15:00", cpus=1, mem=2, nodes=1,
                        qos="regular"):

        decon_qtl_path = "/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/2020-11-20-decon-QTL/Decon-eQTL-v1.4-jar-with-dependencies.jar"

        lines = ["#!/bin/bash",
                 "#SBATCH --job-name={}".format(base_name + job_name),
                 "#SBATCH --output={}".format(os.path.join(jobs_outdir, base_name + job_name + ".out")),
                 "#SBATCH --error={}".format(os.path.join(jobs_outdir, base_name + job_name + ".out")),
                 "#SBATCH --time={}".format(time),
                 "#SBATCH --cpus-per-task {}".format(cpus),
                 "#SBATCH --mem {}gb".format(mem),
                 "#SBATCH --nodes {}".format(nodes),
                 "#SBATCH --qos={}".format(qos),
                 "",
                 "module load Java",
                 "",
                 "java -jar -XX:ParallelGCThreads={} -Xms{}M -Xmx{}M {} \\".format(ParallelGCThreads, initial_heap_size, max_heap_size, decon_qtl_path),
                 "     --outfolder {}/ \\".format(results_outdir),
                 "     --cellcount {} \\".format(cellcount),
                 "     --expression {} \\".format(expression),
                 "     --snpsToTest {} \\".format(snps_to_test),
                 "     --genotype {} \\".format(genotype),
                 "     -gc one",
                 ""]

        jobfile_path = os.path.join(jobsdir, base_name + job_name + ".sh")
        with open(jobfile_path, "w") as f:
            for line in lines:
                f.write(line + "\n")
        f.close()
        print("\tSaved jobfile: {}".format(os.path.basename(jobfile_path)))

    @staticmethod
    def create_perm_orders(n_permutations, cohort_df):
        n_samples = cohort_df.shape[1]
        valid_order = set([x for x in range(n_samples)])
        perm_order_m = np.empty((n_permutations, n_samples), dtype=np.uint8)
        for i in range(perm_order_m.shape[0]):
            sample_indices = np.array([x for x in range(n_samples)])
            for j in range(cohort_df.shape[0]):
                mask = cohort_df.iloc[j, :] == 1
                if np.sum(mask) <= 1:
                    continue
                copy = sample_indices[mask]
                random.shuffle(copy)
                sample_indices[mask] = copy

            if len(set(sample_indices).symmetric_difference(valid_order)) != 0:
                print("Unvalid permutation order.")
                exit()

            perm_order_m[i, :] = sample_indices
        return perm_order_m

    def print_arguments(self):
        print("Arguments:")
        print("  > Deconvolution path: {}".format(self.decon_path))
        print("  > eQTL path: {}".format(self.eqtl_path))
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Cell counts path: {}".format(self.cc_path))
        print("  > Cohort path: {}".format(self.cohort_path))
        print("  > Alpha: {}".format(self.alpha))
        print("  > N permutations: {}".format(self.n_permutations))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
