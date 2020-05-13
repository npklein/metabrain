"""
File:         create_matrices.py
Created:      2020/03/12
Last Changed: 2020/05/13
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
import pickle
import gzip
import os

# Third party imports.

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe
from general.objects.eqtl import Eqtl
from general.objects.group import Group


class CreateMatrices:
    def __init__(self, settings, gte_df, sample_dict, sample_order, eqtl_df,
                 force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param gte_df: DataFrame, the combined GTE files in a dataframe.
        :param sample_dict: dictionary, a dictionary for translating unmasked
                            sampels to the same format.
        :param sample_order: list, order of samples.
        :param eqtl_df: DataFrame, the combined eQTL probe files in a dataframe.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.geno_file = settings["genotype_datafile"]
        self.expr_file = settings["expression_datafile"]
        self.gte_df = gte_df
        self.sample_dict = sample_dict
        self.sample_order = sample_order
        self.eqtl_df = eqtl_df
        self.force = force

        # Prepare an output directories.
        self.outdir = os.path.join(outdir, 'create_matrices')
        prepare_output_dir(self.outdir)

        # Construct the output paths.
        self.geno_outpath = os.path.join(self.outdir, "genotype_table.txt.gz")
        self.alleles_outpath = os.path.join(self.outdir, "genotype_alleles.txt.gz")
        self.expr_outpath = os.path.join(self.outdir, "expression_table.txt.gz")
        self.group_outpath = os.path.join(self.outdir, "groups.pkl")

        # Create empty variable.
        self.complete_expr_matrix = None

    def start(self):
        print("Starting creating matrices.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.geno_outpath) and \
                check_file_exists(self.alleles_outpath) and \
                check_file_exists(self.expr_outpath) and \
                check_file_exists(self.group_outpath) and \
                not self.force:
            print("Skipping step.")
            return

        # Remove the output files.
        for outfile in [self.geno_outpath, self.alleles_outpath,
                        self.expr_outpath, self.group_outpath]:
            if os.path.isfile(outfile):
                print("Removing file: {}.".format(outfile))
                os.remove(outfile)

        # Load the genotype matrix file.
        print("Loading genotype matrix.")
        geno_df = load_dataframe(self.geno_file, header=0, index_col=0)
        allele_df = geno_df.loc[:, ["Alleles", "MinorAllele"]].copy()
        geno_df = geno_df.rename(columns=self.sample_dict)
        geno_df = geno_df[self.sample_order]

        # Load the expression matrix file.
        print("Loading expression matrix.")
        expr_df = load_dataframe(self.expr_file, header=0, index_col=0)
        expr_df = expr_df.rename(columns=self.sample_dict)
        self.complete_expr_matrix = expr_df[self.sample_order]

        # Construct the genotype / expression matrices.
        print("Constructing matrices.")
        geno_str_buffer = ["-" + "\t" + "\t".join(self.sample_order) + "\n"]
        expr_str_buffer = ["-" + "\t" + "\t".join(self.sample_order) + "\n"]
        allele_str_buffer = ["-" + "\t" + "\t".join(list(allele_df.columns)) + "\n"]

        saved_profile_genes = []
        groups = []
        new_group_id = 0
        n_snps = self.eqtl_df.shape[0]
        for i, row in self.eqtl_df.iterrows():
            if (i % 250 == 0) or (i == (n_snps - 1)):
                print("\tProcessing {}/{} "
                      "[{:.2f}%]".format(i,
                                         (n_snps - 1),
                                         (100 / (n_snps - 1)) * i))

                # Write output files.
                self.write_buffer(self.geno_outpath, geno_str_buffer)
                geno_str_buffer = []

                self.write_buffer(self.expr_outpath, expr_str_buffer)
                expr_str_buffer = []

                self.write_buffer(self.alleles_outpath, allele_str_buffer)
                allele_str_buffer = []

            # Get the row info.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]

            # Used for development.
            # snp_name = "10:100145864:rs4919426:T_C"
            # probe_name = "ENSG00000000003.15"
            # End used for development.

            # Get the genotype.
            genotype = geno_df.loc[[snp_name], :]
            if (len(genotype.index)) != 1:
                print("SNP: {} gives 0 or >1 "
                      "genotypes.".format(snp_name))
                continue
            geno_str = snp_name + "\t" + "\t".join(
                genotype.iloc[0, :].astype(str).values) + "\n"
            geno_str_buffer.append(geno_str)

            # Get the alleles.
            alleles = allele_df.loc[[snp_name], :]
            if (len(alleles.index)) != 1:
                print("SNP: {} gives 0 or >1 "
                      "alleles.".format(snp_name))
                continue
            allele_str = "{}\t{}\t{}\n".format(snp_name,
                                               alleles.iloc[0]["Alleles"],
                                               alleles.iloc[0][
                                                   "MinorAllele"])
            allele_str_buffer.append(allele_str)

            # Get the expression.
            expression = self.complete_expr_matrix.loc[[probe_name], :]
            if (len(expression.index)) != 1:
                print("Probe: {} gives 0 or >1 expression "
                      "profiles.".format(probe_name))
                continue
            expr_str = probe_name + "\t" + "\t".join(
                expression.iloc[0, :].astype(str).values) + "\n"
            expr_str_buffer.append(expr_str)

            # Create an eQTL object.
            new_eqtl = Eqtl(snp_name, i, genotype, expression)

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
                new_group = Group(new_group_id, samples)
                new_group.add_eqtl(new_eqtl)
                groups.append(new_group)
                new_group_id = new_group_id + 1

        # Write output files.
        if geno_str_buffer:
            self.write_buffer(self.geno_outpath, geno_str_buffer)

        if expr_str_buffer:
            self.write_buffer(self.expr_outpath, expr_str_buffer)

        if allele_str_buffer:
            self.write_buffer(self.alleles_outpath, allele_str_buffer)

        # Pickle the groups.
        print("Writing group pickle file.")
        with open(self.group_outpath, "wb") as f:
            pickle.dump(groups, f)

        # Remove old dataframes.
        del geno_df, expr_df

    @staticmethod
    def write_buffer(filename, buffer):
        """
        Method for writing a list of strings to a gzipped file,

        :param filename: string, the output file path.
        :param buffer: list, the lines of strings to write.
        """
        # Write output files.
        if os.path.isfile(filename):
            mode = 'ab'
        else:
            mode = 'wb'

        with gzip.open(filename, mode) as f:
            for line in buffer:
                f.write(line.encode())
        f.close()

    def clear_variables(self):
        self.geno_file = None
        self.gte_df = None
        self.sample_dict = None
        self.sample_order = None
        self.eqtl_df = None
        self.force = None

    def get_geno_outpath(self):
        return self.geno_outpath

    def get_alleles_outpath(self):
        return self.alleles_outpath

    def get_expr_file(self):
        return self.expr_file

    def get_expr_outpath(self):
        return self.expr_outpath

    def get_complete_expr_matrix(self):
        return self.complete_expr_matrix

    def get_group_outpath(self):
        return self.group_outpath

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype input file: {}".format(self.geno_file))
        print("  > Expression input file: {}".format(self.expr_file))
        print("  > Genotype output path: {}".format(self.geno_outpath))
        print("  > Alleles output path: {}".format(self.alleles_outpath))
        print("  > Expression output path: {}".format(self.expr_outpath))
        print("  > Groups pickle output path: {}".format(self.group_outpath))
        print("  > GTE input shape: {}".format(self.gte_df.shape))
        print("  > eQTL input shape: {}".format(self.eqtl_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("  > Force: {}".format(self.force))
        print("")
