"""
File:         covariates_explained_by_others.py
Created:      2020/04/15
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
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from colour import Color
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.
from general.utilities import prepare_output_dir


class CovariatesExplainedByOthers:
    def __init__(self, dataset, outdir):
        """
        The initializer for the class.

        :param dataset: Dataset, the input data.
        :param outdir: string, the output directory.
        """
        self.outdir = os.path.join(outdir, 'covariates_explained_by_others')
        prepare_output_dir(self.outdir)

        # Extract the required data.
        print("Loading data")
        self.cov_df = dataset.get_cov_df()
        self.colormap = self.create_color_map()

    @staticmethod
    def create_color_map():
        """
        """
        palette = list(Color("#8ABBDB").range_to(Color("#344A5A"), 101))
        colors = [str(x).upper() for x in palette]
        values = [x / 100 for x in list(range(101))]
        color_map = {}
        for val, col in zip(values, colors):
            color_map[val] = col
        return color_map

    def start(self):
        print("Plotting if covariates are explained by other covariates.")
        self.print_arguments()
        r2_df = self.model(self.cov_df)

        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None):
            print(r2_df)

        self.plot(r2_df, self.outdir)

    def model(self, cov_df):
        print("Modelling each covariate with a linear model of the others.")
        r2_data = []
        for index in cov_df.index:
            X = cov_df.loc[cov_df.index != index, :].T
            y = cov_df.loc[index, :]
            score = self.create_model(X, y)
            color = self.colormap[round(score, 2)]
            r2_data.append([index, score, color])

        r2_df = pd.DataFrame(r2_data, columns=["index", "value", "color"])

        return r2_df

    @staticmethod
    def create_model(X, y):
        """
        Method for creating a multilinear model.

        :param X: DataFrame, the matrix with rows as samples and columns as
                             dimensions.
        :param y: Series, the outcome values.
        :return degrees_freedom: int, the degrees of freedom of this model.
        :return residual_squared_sum: float, the residual sum of squares of this
                                      fit.
        """
        # Create the model.
        regressor = LinearRegression()
        regressor.fit(X, y)
        y_hat = regressor.predict(X)

        # Calculate the statistics of the model.
        score = r2_score(y, y_hat)

        return score

    @staticmethod
    def plot(df, outdir):
        print("Plotting")

        indices = [(0, 18, "Tech. Cov.", ""),
                   (18, 22, "MDS", "MDS"),
                   (22, 40, "Cohorts", ""),
                   (40, 41, "Sex", ""),
                   (41, 91, "PCs", "Comp"),
                   (91, 93, "PCA", ""),
                   (93, 118, "McKenzie\nMG", "McKenzie_"),
                   (118, 123, "CellMap\nPCA", "CellMapPCA_"),
                   (123, 128, "CellMap\nNMF", "CellMapNMF_"),
                   (128, 133, "CellMap\nNNLS", "CellMapNNLS_")]

        gridspec_kw = {"height_ratios": [x[1] - x[0] for x in indices],
                       "width_ratios": [0.3, 0.7]}

        sns.set(style="ticks", color_codes=True)
        fig, axes = plt.subplots(ncols=2, nrows=len(indices),
                                 figsize=(9, 28), gridspec_kw=gridspec_kw)
        plt.subplots_adjust(top=0.95, bottom=0.05, wspace=0.1, hspace=0.2)

        for i in range(len(indices)):
            print("\tPlotting axes[{}, 1]".format(i))
            axes[i, 0].set_axis_off()
            ax = axes[i, 1]
            sns.despine(fig=fig, ax=ax)
            (a, b, ylabel, remove) = indices[i]
            xlabel = ""
            if i == (len(indices) - 1):
                xlabel = "MLR R2"

            subset = df.iloc[a:b, :]

            for i in range(0, 105, 5):
                alpha = 0.025
                if i % 10 == 0:
                    alpha = 0.15
                ax.axvline(i / 100, ls='-', color="#000000",
                           alpha=alpha, zorder=-1)

            sns.barplot(x="value", y="index", data=df.iloc[a:b, :],
                        palette=subset["color"], orient="h", ax=ax)

            new_ylabels = [x.replace(remove, '') for x in subset["index"]]
            ax.set_yticklabels(new_ylabels, fontsize=10)
            ax.set_ylabel(ylabel, fontsize=16, fontweight='bold')
            ax.set_xlabel(xlabel, fontsize=16, fontweight='bold')
            ax.set(xlim=(0, 1))

        fig.align_ylabels(axes[:, 1])
        fig.suptitle('Variance Explained by other Covariates', fontsize=25, fontweight='bold')
        fig.savefig(os.path.join(outdir, "covariates_explained_by_others.png"))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Covariate matrix shape: {}".format(self.cov_df.shape))
        print("  > Output directory: {}".format(self.outdir))
        print("")
