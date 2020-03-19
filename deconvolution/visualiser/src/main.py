"""
File:         main.py
Created:      2020/03/13
Last Changed: 2020/03/18
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
from general.utilities import prepare_output_dir
from general.local_settings import LocalSettings
from general.objects.dataset import Dataset
from .figures.simple_eqtl_effect import SimpleeQTLEffect
from .figures.inter_zscore_bars import InterZscoreBars
from .figures.inter_zscore_dist import InterZscoreDist
from .figures.inter_zscore_clustermap import InterZscoreClusterMap
from .figures.inter_eqtl_zscore_bars import IntereQTLZscoreBars
from .figures.inter_zscore_marker_genes import InterZscoreMarkerGenes
from .figures.inter_eqtl_effect import IntereQTLEffect
from .figures.inter_eqtl_effect_marker_genes import IntereQTLEffectMarkerGenes
from .figures.inter_eqtl_effect_marker_vs_comp import IntereQTLEffectMarkerVSComp


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
        # Define the current directory.
        current_dir = str(Path(__file__).parent.parent)

        # Load the LocalSettings singelton class.
        self.settings = LocalSettings(current_dir, settings_file)

        # Load the variables.
        self.plots = plots

        # Prepare an output directory.
        self.outdir = os.path.join(current_dir,
                                   self.settings.get_setting('output_dir'),
                                   self.settings.get_setting('group_dir'))
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
        if ('simple_eqtl_effect' in self.plots) or ('all' in self.plots):
            print("\n### SIMPLE EQTL EFFECT ###\n")
            sef = SimpleeQTLEffect(dataset=ds,
                                   outdir=self.outdir)
            sef.start()
            del sef

        if ('inter_zscore_bars' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION Z-SCORE BARPLOT ###\n")
            izb = InterZscoreBars(dataset=ds,
                                  outdir=self.outdir)
            izb.start()
            del izb

        if ('inter_zscore_dist' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION Z-SCORE DISTRIBUTION PLOT ###\n")
            izd = InterZscoreDist(dataset=ds,
                                  outdir=self.outdir)
            izd.start()
            del izd

        if ('inter_zscore_marker_genes' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION Z-SCORE CLUSTERMAP ###\n")
            izmg = InterZscoreMarkerGenes(dataset=ds,
                                          outdir=self.outdir)
            izmg.start()
            del izmg

        if ('inter_zscore_clustermap' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION Z-SCORE CLUSTERMAP ###\n")
            izcp = InterZscoreClusterMap(dataset=ds,
                                         outdir=self.outdir)
            izcp.start()
            del izcp

        if ('inter_eqtl_zscore_bars' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL Z-SCORE BARS ###\n")
            iezb = IntereQTLZscoreBars(dataset=ds,
                                       outdir=self.outdir)
            iezb.start()
            del iezb

        if ('inter_eqtl_effect' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT ###\n")
            iee = IntereQTLEffect(dataset=ds,
                                  outdir=self.outdir)
            iee.start()
            del iee

        if ('inter_eqtl_effect_marker_genes' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT MARKER GENES ###\n")
            ieemg = IntereQTLEffectMarkerGenes(dataset=ds,
                                               outdir=self.outdir)
            ieemg.start()
            del ieemg

        if ('inter_eqtl_effect_marker_vs_comp' in self.plots) or ('all' in self.plots):
            print("\n### INTERACTION EQTL EFFECT MARKER GENES VS COMPS ###\n")
            ieemvc = IntereQTLEffectMarkerVSComp(dataset=ds,
                                                 outdir=self.outdir)
            ieemvc.start()
            del ieemvc

    def print_arguments(self):
        print("Arguments:")
        print("  > Output directory: {}".format(self.outdir))
        print("  > Plots: {}".format(self.plots))
        print("")
