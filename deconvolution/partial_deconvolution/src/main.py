"""
File:         main.py
Created:      2020/06/29
Last Changed: 2020/09/03
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
import os

# Third party imports.

# Local application imports.
from .data_loader import DataLoader
from .data_preprocessor import DataPreprocessor
from .perform_deconvolution import PerformDeconvolution
from .data_comparitor import DataComparitor
from .visualiser import Visualiser


class Main:
    def __init__(self, settings, outdir, visualise, plot_ids):
        self.settings = settings
        self.visualise = visualise
        self.plot_ids = plot_ids

        # Prepare output directory.
        current_dir = str(Path(__file__).parent.parent)
        outpath = os.path.join(current_dir, outdir)
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        self.settings.set_output_path(outpath)

    def start(self):
        # Load the data.
        print("### Loading")
        dl = DataLoader(self.settings)
        dl.work()
        dl.print_info()
        self.settings.set_real_info_per_celltype(dl.get_info_per_celltype())

        # Save.
        self.settings.save_settings()

        # Preprocessing.
        print("### Preprocessing")
        dp = DataPreprocessor(settings=self.settings,
                              raw_signature=dl.get_signature(),
                              raw_expression=dl.get_expression(),
                              cohorts=dl.get_cohorts())
        dp.work()
        dp.print_info()
        self.settings.set_filter_shape_diff(dp.get_shape_diff())
        self.settings.set_sign_shift(dp.get_sign_shift())
        self.settings.set_expr_shift(dp.get_expr_shift())

        # Partial deconvolution.
        print("### Deconvoluting")
        pf = PerformDeconvolution(settings=self.settings,
                                  signature=dp.get_signature(),
                                  expression=dp.get_expression())
        pf.work()
        pf.print_info()
        self.settings.set_avg_residuals(pf.get_avg_residuals())
        self.settings.set_pred_info_per_celltype(pf.get_info_per_celltype())

        # Comparison.
        print("### Comparing")
        dc = DataComparitor(settings=self.settings,
                            deconvolution=pf.get_deconvolution(),
                            ground_truth=dl.get_ground_truth())
        dc.work()
        dc.print_info()
        self.settings.set_comparison_n_samples(dc.get_n_samples())
        self.settings.set_comparison_rss(dc.get_rss())

        # Save.
        self.settings.save_settings()

        # Visualising profile.
        if self.visualise:
            print("### Visualising")
            v = Visualiser(settings=self.settings,
                           signature=dp.get_signature(),
                           expression=dp.get_expression(),
                           deconvolution=pf.get_deconvolution(),
                           ground_truth=dl.get_ground_truth(),
                           comparison=dc.get_comparison())
            # v.plot_profile_clustermap()
            # v.plot_profile_stripplot()
            # v.plot_profile_boxplot()
            # v.plot_deconvolution_clustermap()
            # v.plot_deconvolution_per_sample()
            # #v.plot_deconvolution_distribution()
            # v.plot_deconvolution_boxplot()
            # #v.plot_ground_truth_distribution()
            # v.plot_ground_truth_boxplot()
            # v.plot_prediction_comparison()

            for plot_id in self.plot_ids:
                print("Plotting {}".format(plot_id))
                v.plot_violin_comparison(plot_id, dl.get_sample_to_info_dict(plot_id))