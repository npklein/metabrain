"""
File:         cohort_clustermap.py
Created:      2020/06/02
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

# Local application imports.
from general.utilities import prepare_output_dir


class CohortClustermap:
    def __init__(self, dataset, outdir, extension):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        :param extension: str, the output figure file type extension.
        """
        self.outdir = os.path.join(outdir, 'cohort_clustermap')
        prepare_output_dir(self.outdir)
        self.extension = extension

        # Extract the required data.
        print("Loading data")
        self.cov_df = dataset.get_cov_df()
        self.cohorts = ["AMPAD-MSBB-V2-AFR",
                        "CMC_HBCC_set1-AFR",
                        "CMC_HBCC_set2-AFR",
                        "CMC-AFR",
                        "ENA-AFR",
                        "LIBD_1M-AFR",
                        "LIBD_h650-AFR",
                        "ENA-EAS",
                        "AMPAD-MAYO-V2-EUR",
                        "AMPAD-MSBB-V2-EUR",
                        "BrainGVEX-V2-EUR",
                        "CMC_HBCC_set2-EUR",
                        "CMC_HBCC_set3-EUR",
                        "CMC-EUR",
                        "ENA-EUR",
                        "GTEx-EUR",
                        "GVEX-EUR",
                        "LIBD_1M-EUR",
                        "LIBD_h650-EUR",
                        "NABEC-H550-EUR",
                        "NABEC-H610-EUR",
                        "TargetALS-EUR",
                        "UCLA_ASD-EUR",
                        "OTHER"
                        ]

    def start(self):
        print("Plotting convariate comparison.")
        self.print_arguments()
        cohorts = self.cov_df.loc[self.cohorts, :].copy()
        self.plot(cohorts, self.outdir, self.extension)

    @staticmethod
    def plot(df, outdir, extension):
        sns.set(color_codes=True)
        g = sns.clustermap(df, center=0.5, cmap="RdBu_r",
                           vmin=0, vmax=1,
                           yticklabels=True, xticklabels=False,
                           dendrogram_ratio=(.1, .1),
                           figsize=(12, (.2 * (len(df.index)))))
        g.cax.set_visible(False)
        plt.setp(
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(),
                                         fontsize=10))
        g.savefig(os.path.join(outdir, "cohorts_clustermap.{}".format(extension)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Cohorts: {}".format(self.cohorts))
        print("  > Output directory: {}".format(self.outdir))
        print("")
