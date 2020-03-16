"""
File:         inter_eqtl_zscore_bars.py
Created:      2020/03/16
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
import seaborn as sns
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from colour import Color
import scipy.stats as st

# Local application imports.
from src.utilities import prepare_output_dir


class IntereQTLZscoreBars:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'inter_eqtl_zscore_bars')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        self.eqtl_df = dataset.get_eqtl_df()
        self.inter_df = dataset.get_inter_df()

    def start(self):
        print("Plotting interaction eQTL z-score barplots.")
        self.print_arguments()

        # Calculate the z-score cutoff.
        z_score_cutoff = st.norm.ppf(
            0.05 / (self.inter_df.shape[0] * self.inter_df.shape[1]) / 2)

        i = 0
        for index, row in self.eqtl_df.iterrows():
            # Extract the usefull information from the row.
            snp_name = row["SNPName"]
            probe_name = row["ProbeName"]
            hgnc_name = row["HGNCName"]
            eqtl_type = row["CisTrans"]

            # Check if the SNP has an interaction effect.
            interaction_effect = self.inter_df.iloc[:, i].to_frame()
            # if interaction_effect.columns[0] in ["7:144177389:rs6464583:A_G_G", "9:32886201:rs10971181:T_C_C", "22:16717912:rs2078647:A_G_G", "15:33985176:rs11856529:A_G_G"]:

            interaction_effect.reset_index(inplace=True)
            interaction_effect.columns = ["index", "zscore"]
            interaction_effect = interaction_effect.reindex(
                interaction_effect["zscore"].sort_values().index)

            self.plot(i, snp_name, probe_name, hgnc_name, eqtl_type,
                      z_score_cutoff, interaction_effect, self.outdir)

        i += 1

    @staticmethod
    def plot(i, snp_name, probe_name, hgnc_name, eqtl_type, z_score_cutoff,
             df, outdir):
        sns.set(rc={'figure.figsize': (12, (.2 * (len(df.index))))})
        fig, ax = plt.subplots()
        g = sns.barplot(x="zscore", y="index", data=df, palette="coolwarm",
                        orient="h")
        g.text(0.5, 1.025,
               '{} {}-eQTL Interactions'.format(hgnc_name, eqtl_type),
               fontsize=22, weight='bold', ha='center', va='bottom',
               transform=ax.transAxes)
        g.text(0.5, 1.01,
               'z-score cutoff = +/- {:.2f}'.format(abs(z_score_cutoff)),
               fontsize=20, alpha=0.75, ha='center', va='bottom',
               transform=ax.transAxes)
        g.set_ylabel('covariate',
                     fontsize=18,
                     fontweight='bold')
        g.set_xlabel('z-score',
                     fontsize=18,
                     fontweight='bold')
        plt.axvline(x=z_score_cutoff, ls='--', c='red')
        plt.axvline(x=-1 * z_score_cutoff, ls='--', c='red')
        ax.tick_params(labelsize=10)
        ax.set_yticks(range(len(df.index)))
        ax.set_yticklabels(df["index"])
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,
                                 "{}_inter_eqtl_bars_{}_{}_{}.png".format(i,
                                                                          snp_name,
                                                                          probe_name,
                                                                          hgnc_name)))
        plt.close()


    def print_arguments(self):
        print("Arguments:")
        print("  > eQTL matrix shape: {}".format(self.eqtl_df.shape))
        print("  > Interaction matrix shape: {}".format(self.inter_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
