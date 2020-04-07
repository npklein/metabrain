"""
File:         combine_gte_files.py
Created:      2020/04/07
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
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

# Local application imports.
from general.utilities import prepare_output_dir, check_file_exists
from general.df_utilities import load_dataframe, save_dataframe


class PerformCelltypePCA:
    def __init__(self, settings, profile_file, profile_df, ct_expr_file,
                 force, outdir):
        """
        The initializer for the class.


        :param settings: string, the settings.
        :param profile_file: string, the datafile contaioning the celltype
                             profile.
        :param profile_df: DataFrame, the celltype profile.
        :param ct_expr_file: string, the datafile containing expression
                             of the celltype profiles.
        :param force: boolean, whether or not to force the step to redo.
        :param outdir: string, the output directory.
        """
        self.profile_file = profile_file
        self.profile_df = profile_df
        self.ct_expr_file = ct_expr_file
        self.force = force

        # Prepare an output directory.
        self.outdir = os.path.join(outdir, 'perform_celltype_pca')
        prepare_output_dir(self.outdir)
        self.outpath = os.path.join(self.outdir, "celltype_pcs.txt.gz")

        # Create empty variables.
        self.celltype_pcs = None

    def start(self):
        print("Starting PCA of celltype profile expression.")
        self.print_arguments()

        # Check if output file exist.
        if check_file_exists(self.outpath) and not self.force:
            print("Skipping step, loading result.")
            self.celltype_pcs = load_dataframe(inpath=self.outpath,
                                               header=0, index_col=0)
        else:
            self.celltype_pcs = self.perform_pca()
            self.save()

    def perform_pca(self):
        # Load the expression data.
        print("Loading celltype expression data.")
        ct_expr_df = load_dataframe(inpath=self.ct_expr_file,
                                    header=0, index_col=0)

        if self.profile_df is None:
            # Load the celltype profile file.
            print("Loading cell type profile matrix.")
            self.profile_df = load_dataframe(self.profile_file,
                                             header=0, index_col=0)

        # Find the genes specific to each celltype.
        gene_celltypes = self.normalize(self.profile_df).idxmax(axis=1)

        # Construct a dataframe of the first component of each celltype
        # subset expression profile.
        pc_data = []
        print("Performing PCA")
        for celltype in self.profile_df.columns:
            print("\tWorking on: {}.".format(celltype))
            ct_genes = gene_celltypes[gene_celltypes == celltype].index
            ct_expr = ct_expr_df.loc[ct_expr_df.index.isin(ct_genes), :]

            # perform PCA over the expression of these genes.
            component = self.get_first_component(ct_expr)
            component_values = [x[0] for x in list(component)]
            pc_data.append(component_values)

        # Create a dataframe.
        celltype_pcs = pd.DataFrame(pc_data,
                                    index=["{}_PC1".format(x) for x in self.profile_df.columns],
                                    columns=ct_expr_df.columns)

        return celltype_pcs

    def save(self):
        save_dataframe(df=self.celltype_pcs, outpath=self.outpath,
                       index=True, header=True)

    @staticmethod
    def normalize(df):
        df = df.loc[df.std(axis=1) > 0, :]
        df = df.T
        df = (df-df.mean()) / df.std()
        df = df.T
        return df

    @staticmethod
    def get_first_component(X):
        corr_matrix = np.dot(X.T, X) / (X.shape[0] - 1)

        pca = PCA(n_components=1)
        pca.fit(corr_matrix)
        print("\t  Explained variance ratio: {}".format(pca.explained_variance_ratio_[0]))
        print("\t  Singular values: {}".format(pca.singular_values_[0]))
        return pca.transform(corr_matrix)

    def get_celltype_pcs(self):
        return self.celltype_pcs

    def clear_variables(self):
        self.profile_file = None
        self.profile = None
        self.ct_expr_file = None
        self.force = None

    def print_arguments(self):
        print("Arguments:")
        if self.profile_df is not None:
            print("  > Celltype profile: {}".format(self.profile_df.shape))
        else:
            print("  > Celltype profile input file: {}".format(self.profile_file))
        print("  > Celltype expression input path: {}".format(self.ct_expr_file))
        print("  > Celltype PCAs outpath file: {}".format(self.outpath))
        print("  > Force: {}".format(self.force))
        print("")
