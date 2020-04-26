"""
File:         inter_eqtl_effect_celltype.py
Created:      2020/04/20
Last Changed: 2020/04/26
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
import matplotlib
matplotlib.use('Agg')
from venn import venn
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class IntereQTLEffectCelltype:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_effect_celltype')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.eqtl_ct_df = dataset.get_eqtl_ct_df()
        self.celltypes = dataset.get_celltypes()
        self.cellmap_methods = dataset.get_cellmap_methods()
        self.marker_genes = dataset.get_marker_genes()
        self.colormap = dataset.get_colormap()

    def start(self):
        print("Plotting cell type mediated interaction eQTLs")
        self.print_arguments()

        celltype_mediated_eqtls = []

        print("Plotting marker genes.")
        for celltype in self.celltypes:
            df = self.eqtl_ct_df.loc[:, self.eqtl_ct_df.columns.str.startswith(self.marker_genes + celltype + "_")]
            df.index = self.eqtl_ct_df["SNPName"]
            data = {}
            for column in df:
                gene_name = column.split("_")[2]
                eqtls = set(df.loc[df[column] == 1, column].index.values)
                celltype_mediated_eqtls.extend(eqtls)
                data[gene_name] = eqtls
            self.plot_venn(data,
                           "{}{}".format(self.marker_genes, celltype),
                           self.outdir)

        print("Plotting deconvolution methods.")
        for celltype in self.celltypes:
            data = {}

            # Add the marker genes.
            marker_df = self.eqtl_ct_df.loc[:, self.marker_genes + celltype]
            marker_df.index = self.eqtl_ct_df["SNPName"]
            eqtls = set(marker_df.loc[marker_df == 1].index.values)
            celltype_mediated_eqtls.extend(eqtls)
            data[self.marker_genes.replace("_", "")] = eqtls

            method_celltype = celltype
            title = celltype
            if method_celltype == "Microglia":
                method_celltype = "Macrophage"
                title = method_celltype + "_" + celltype

            # Add the deconvolution methods.
            for (prefix, suffix) in self.cellmap_methods:
                df = self.eqtl_ct_df.loc[:, prefix + method_celltype + suffix]
                df.index = self.eqtl_ct_df["SNPName"]
                eqtls = set(df.loc[df == 1].index.values)
                celltype_mediated_eqtls.extend(eqtls)
                data[prefix.replace("_", "") + suffix] = eqtls

            # Plot.
            self.plot_venn(data,
                           title,
                           self.outdir)

        print("Plotting all eQTLs.")
        data = {"cis-eQTL" : set(self.eqtl_ct_df["SNPName"].values),
                "cell type mediated": set(celltype_mediated_eqtls)}
        self.plot_venn(data,
                       "all_cis-eQTLs",
                       self.outdir)

    @staticmethod
    def plot_venn(data, title, outdir):
        fig, ax = plt.subplots(1, 1, figsize=(12, 9))
        venn(data, fontsize=14, legend_loc="upper right", ax=ax)
        ax.text(0.5, 1.05,
                title,
                fontsize=20, weight='bold', ha='center', va='bottom',
                transform=ax.transAxes)
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}.png".format(title)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_ct_df.shape))
        print("  > Celltypes: {}".format(self.celltypes))
        print("  > CellMap Methods: {}".format(self.cellmap_methods))
        print("  > Marker Genes: {}".format(self.marker_genes))
        print("  > Output directory: {}".format(self.outdir))
        print("")
