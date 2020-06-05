"""
File:         combine_eqtlprobes.py
Created:      2020/03/12
Last Changed: 2020/06/03
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
import pandas as pd

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe, save_dataframe


class CombineEQTLProbes:
    def __init__(self, settings, disease, force, outdir):
        """
        The initializer for the class.

        :param settings: string, the settings.
        :param disease: string, the name of the disease to analyse.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.indir = settings["input_directory"]
        self.iter_dirname = settings["iteration_dirname"]
        self.n_iterations = settings["iterations"]
        self.in_filename = settings["in_filename"]
        self.snp_to_gwasid_filename = settings["snp_to_gwasid_filename"]
        self.gwasid_to_trait_filename = settings["gwasid_to_trait_filename"]
        self.disease = disease
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'combine_eqtlprobes')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "eQTLprobes_combined.txt.gz")

        # Declare variables.
        self.eqtl_probes = None

    def start(self):
        print("Starting combining eQTL probe files.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.eqtl_probes = load_dataframe(inpath=self.outpath, header=0,
                                              index_col=False)
        else:
            # Load each GTE file.
            print("Loading eQTLprobes files.")
            combined_eqtl_probes = self.combine_files()
            if self.disease != "" and self.disease is not None:
                print("Filtering on trait: {}".format(self.disease))
                combined_eqtl_probes = self.filter_on_trait(combined_eqtl_probes)
            self.eqtl_probes = combined_eqtl_probes
            self.save()

    def combine_files(self):
        combined = None
        for i in range(1, self.n_iterations+1):
            infile = os.path.join(self.indir, self.iter_dirname + str(i),
                                  self.in_filename)
            df = load_dataframe(inpath=infile, header=0, index_col=False)
            df["Iteration"] = i
            if combined is None:
                combined = df
            else:
                combined = pd.concat([combined, df], axis=0, ignore_index=True)

        # Remove duplicate entries.
        combined.drop_duplicates(inplace=True)

        return combined

    def filter_on_trait(self, df):
        tmp1 = load_dataframe(inpath=self.gwasid_to_trait_filename, header=0,
                              index_col=False)
        gwas_to_trait = pd.Series(tmp1["Trait"].values, index=tmp1["ID"]).to_dict()
        del tmp1

        gwas_map = {}
        disease_map = {}
        tmp2 = load_dataframe(inpath=self.snp_to_gwasid_filename, header=0,
                              index_col=False, low_memory=False)

        for index, row in tmp2.iterrows():
            rs = row["RsID"]
            id = row["ID"]

            gwasses = gwas_map.get(rs)
            if gwasses is None:
                gwasses = id
            else:
                gwasses = "{}, {}".format(gwasses, id)
            gwas_map[rs] = gwasses

            diseases = disease_map.get(rs)
            if id in gwas_to_trait.keys():
                trait = gwas_to_trait.get(id)
                if diseases is None:
                    diseases = trait
                else:
                    diseases = "{}, {}".format(diseases, trait)
            disease_map[rs] = diseases

        df["GWASIDS"] = df["SNPName"].map(gwas_map, na_action="")
        df["Trait"] = df["SNPName"].map(disease_map, na_action="")

        # Subset.
        df.dropna(subset=['Trait'], inplace=True)
        df = df[df['Trait'].str.contains(self.disease, case=False)]
        df.reset_index(drop=True, inplace=True)

        return df

    def save(self):
        save_dataframe(df=self.eqtl_probes, outpath=self.outpath,
                       index=False, header=True)

    def clear_variables(self):
        self.indir = None
        self.iter_dirname = None
        self.n_iterations = None
        self.in_filename = None
        self.force = None

    def get_outpath(self):
        return self.outpath

    def get_eqtlprobes(self):
        return self.eqtl_probes

    def print_arguments(self):
        print("Arguments:")
        print("  > Input directory: {}".format(self.indir))
        print("  > Iteration directory: {}".format(self.iter_dirname))
        print("  > N. Iterations: {}".format(self.n_iterations))
        print("  > Input filename: {}".format(self.in_filename))
        print("  > Output path: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
