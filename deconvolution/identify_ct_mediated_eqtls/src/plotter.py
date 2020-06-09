"""
File:         plotter.py
Created:      2020/06/09
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
import matplotlib
import itertools
matplotlib.use('Agg')
import upsetplot as up
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class Plotter:
    def __init__(self, df, exclude, total, outdir, extension="png"):
        self.df = df.loc[~df["Covariate"].isin(exclude), :]
        self.total = total
        self.outdir = os.path.join(outdir, 'plots')
        self.extension = extension

        prepare_output_dir(self.outdir)

        # Set the right pdf font for exporting.
        if self.extension == "pdf":
            matplotlib.rcParams['pdf.fonttype'] = 42

    def plot(self):
        print("Plotting interaction upsetplot.")
        data = {}
        for covariate in self.df["Covariate"].unique():
            subset = self.df.loc[self.df["Covariate"] == covariate, :]
            data[covariate] = set(subset["Index"])

        self.upsetplot(data, self.outdir, self.extension)
        self.plot_pie(self.total, self.df.shape[0], self.outdir, self.extension)

    def upsetplot(self, data, outdir, extension):
        counts = self.count(data)
        up.plot(counts, sort_by='cardinality', show_counts=True)
        plt.savefig(os.path.join(outdir, "upsetplot.{}".format(extension)))
        plt.close()

    @staticmethod
    def count(input_data):
        combinations = []
        cols = list(input_data.keys())
        for i in range(1, len(cols) + 1):
            combinations.extend(list(itertools.combinations(cols, i)))

        abbreviations = {"neuron": "neuro",
                         "oligodendrocyte": "oligo",
                         "endothelialcell": "endo",
                         "microglia": "micro",
                         "macrophage": "macro",
                         "astrocyte": "astro"}

        abbr_cols = []
        for col in cols:
            if col in abbreviations.keys():
                abbr_cols.append(abbreviations[col])
            else:
                abbr_cols.append(col)

        indices = []
        data = []
        for combination in combinations:
            index = []
            for col in cols:
                if col in combination:
                    index.append(True)
                else:
                    index.append(False)

            background = set()
            for key in cols:
                if key not in combination:
                    work_set = input_data[key].copy()
                    background.update(work_set)

            overlap = None
            for key in combination:
                work_set = input_data[key].copy()
                if overlap is None:
                    overlap = work_set
                else:
                    overlap = overlap.intersection(work_set)

            duplicate_set = overlap.intersection(background)
            length = len(overlap) - len(duplicate_set)

            indices.append(index)
            data.append(length)

        s = pd.Series(data,
                      index=pd.MultiIndex.from_tuples(indices, names=abbr_cols))
        s.name = "value"
        return s

    @staticmethod
    def plot_pie(total, part, outdir, extension):
        labels = ['YES', 'NO']
        sizes = [part, total - part]
        explode = (0.1, 0)

        fig, ax = plt.subplots()
        _, _, autotexts = ax.pie(sizes,
                                 explode=explode,
                                 labels=labels,
                                 colors=["#6495ED", "#808080"],
                                 startangle=90,
                                 autopct=lambda pct: "{:.1f}%\n(n={:d})".format(
                                     pct, int(pct / 100. * sum(sizes))),
                                 textprops=dict(fontweight="bold",
                                                fontsize=20))
        for autotext in autotexts:
            autotext.set_color('white')
        ax.axis('equal')

        fig.suptitle("Cell-Type Mediated eQTLs", fontsize=20, fontweight='bold')
        fig.savefig(os.path.join(outdir, "all_cis-eQTLs.{}".format(extension)))
        plt.show()
