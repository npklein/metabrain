"""
File:         main.py
Created:      2020/03/13
Last Changed: 2020/03/17
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
import os

# Third party imports.

# Local application imports.
from src.utilities import get_project_root_dir, prepare_output_dir
from src.local_settings import LocalSettings
from src.objects.dataset import Dataset
from src.figures.simple_eqtl_effect import SimpleeQTLEffect
from src.figures.inter_zscore_bars import InterZscoreBars
from src.figures.inter_zscore_dist import InterZscoreDist
from src.figures.inter_zscore_clustermap import InterZscoreClusterMap
from src.figures.inter_eqtl_zscore_bars import IntereQTLZscoreBars
from src.figures.inter_zscore_marker_genes import InterZscoreMarkerGenes
from src.figures.inter_eqtl_effect import IntereQTLEffect
from src.figures.inter_eqtl_effect_marker_genes import IntereQTLEffectMarkerGenes


class Main:
    """
    Main: this class is the main class that calls all other functionality.
    """

    def __init__(self, settings_file, plots):
        """
        Initializer of the class.

        :param settings_file: string, the name of the settings file.
        :param plots: list, the names of the plots to create.
        """
        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(settings_file)

        # Load the variables.
        self.plots = plots

        # Prepare an output directory.
        self.outdir = os.path.join(get_project_root_dir(),
                                   self.settings.get_setting('output_dir'))
        prepare_output_dir(self.outdir)

    def start(self):
        """
        The method that serves as the pipeline of the whole program.
        """
        print("Starting visualiser.")
        self.print_arguments()

        # Create the dataset object.
        ds = Dataset(settings=self.settings)

        # Figure 1: a simple eQTL effect.
        print("\n### SIMPLE EQTL EFFECT ###\n")
        if ('simple_eqtl_effect' in self.plots) or ('all' in self.plots):
            sef = SimpleeQTLEffect(dataset=ds,
                                   outdir=self.outdir)
            sef.start()
            del sef

        print("\n### INTERACTION Z-SCORE BARPLOT ###\n")
        if ('inter_zscore_bars' in self.plots) or ('all' in self.plots):
            izb = InterZscoreBars(dataset=ds,
                                  outdir=self.outdir)
            izb.start()
            del izb

        print("\n### INTERACTION Z-SCORE DISTRIBUTION PLOT ###\n")
        if ('inter_zscore_dist' in self.plots) or ('all' in self.plots):
            izd = InterZscoreDist(dataset=ds,
                                  outdir=self.outdir)
            izd.start()
            del izd

        print("\n### INTERACTION Z-SCORE CLUSTERMAP ###\n")
        if ('inter_zscore_marker_genes' in self.plots) or ('all' in self.plots):
            izmg = InterZscoreMarkerGenes(dataset=ds,
                                          outdir=self.outdir)
            izmg.start()
            del izmg

        print("\n### INTERACTION Z-SCORE CLUSTERMAP ###\n")
        if ('inter_zscore_clustermap' in self.plots) or ('all' in self.plots):
            izcp = InterZscoreClusterMap(dataset=ds,
                                         outdir=self.outdir)
            izcp.start()
            del izcp

        print("\n### INTERACTION EQTL Z-SCORE BARS ###\n")
        if ('inter_eqtl_zscore_bars' in self.plots) or ('all' in self.plots):
            iezb = IntereQTLZscoreBars(dataset=ds,
                                       outdir=self.outdir)
            iezb.start()
            del iezb

        print("\n### INTERACTION EQTL EFFECT ###\n")
        if ('inter_eqtl_effect' in self.plots) or ('all' in self.plots):
            iee = IntereQTLEffect(dataset=ds,
                                  outdir=self.outdir)
            iee.start()
            del iee

        print("\n### INTERACTION EQTL EFFECT MARKER GENES ###\n")
        if ('inter_eqtl_effect_marker_genes' in self.plots) or ('all' in self.plots):
            ieemg = IntereQTLEffectMarkerGenes(dataset=ds,
                                               outdir=self.outdir)
            ieemg.start()
            del ieemg

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > Plots: {}".format(self.plots))
        print("")
