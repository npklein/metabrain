#!/usr/bin/env python3

"""
File:         pre_process_TPM_expression_matrix.py
Created:      2021/11/23
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
from __future__ import print_function
from pathlib import Path
import argparse
import glob
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Local application imports.

# Metadata
__program__ = "Pre-process TPM Expression Matrix"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax: 
# Calculon
./pre_process_TPM_expression_matrix.py -d /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-01-31-expression-tables/2020-01-31-raw-count-tables/ -g /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/gencode.v32.primary_assembly.annotation-genelengths.txt.gz -s /groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-10-12-deconvolution/deconvolution/data/MetaBrain.allCohorts.2020-02-16.TMM.freeze2dot1.Samples.txt.gz -o MetaBrain_NoENA


# Datasets
TargetALS = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/TargetALS/pipelines/no_patch_chromosomes/results/TargetALS.geneCounts.2019-12-13.txt
Braineac = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/ucl-upload-biogen/pipelines/no_patch_chromosomes/results/Braineac.geneCounts.txt.gz
GTEx = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/GTEx/pipelines/no_patch_chromosomes/results/GTEx.geneCounts.txt 
NABEC = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/NABEC/pipelines/no_patch_chromosomes/results/NABEC.geneCounts.txt.gz 
ENA
BrainGVEx = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/BrainGVEx/results/BrainGVEx.geneCounts.2019-12-16.txt 
BipSeq = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/BipSeq/results/BipSeq.geneCounts.2019-12-17.txt 
UCLA_ASD = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/UCLA_ASD/results/UCLA_ASD.geneCounts.2019-12-15.txt 
CMC_HBCC = Calculon:/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/pipelines/no_patch_chromosomes/CMC_HBCC/results/CMC_HBCC.geneCounts.2020-01-02.txt
CMC = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/CMC/pipelines/no_patch_chromosomes/results/CMC.geneCounts.txt.gz 
MSBB = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/MSBB/results/MSBB.geneCounts.2020-01-02.txt.gz
ROSMAP = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/ROSMAP/results/ROSMAP.geneCounts.2020-01-02.txt.gz  
Mayo_CBE = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/MayoCBE/results/MayoCBE.geneCounts.2020-01-06.txt.gz 
Mayo_TCX = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/AMP_AD/pipelines/no_patch_chromosomes/MayoTCX/results/MayoTCX.geneCounts.2020-01-06.txt.gz 
Brainseq = Gearshift:/groups/umcg-biogen/tmp01/input/rawdata/Brainseq/pipelines/no_patch_chromosomes/results/Brainseq.geneCounts.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')
        self.gene_info_path = getattr(arguments, 'gene_info')
        self.samples_path = getattr(arguments, 'samples')
        outfolder = getattr(arguments, 'outfolder')

        # Set variables.
        outdir = str(Path(__file__).parent.parent)
        self.plot_outdir = os.path.join(outdir, 'pre_process_TPM_expression_matrix', outfolder, 'plot')
        self.file_outdir = os.path.join(outdir, 'pre_process_TPM_expression_matrix', outfolder, 'data')
        for outdir in [self.plot_outdir, self.file_outdir]:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

        self.palette = {
            "TargetALS": "#d5c77a",
            "Braineac": "#808080",
            "GTEx": "#0fa67d",
            "NABEC": "#6d743a",
            "ENA": "#d46727",
            "Brain GVEx": "#48b2e5",
            "BipSeq": "#000000",
            "UCLA ASD": "#f36d2a",
            "CMC HBCC": "#0877b4",
            "CMC": "#eae453",
            "MSBB": "#5cc5bf",
            "ROSMAP": "#6950a1",
            "Mayo_CBE": "#9c9fa0",
            "Mayo_TCX": "#9c9fa0",
            "Brainseq": "#e49d26"
        }

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data matrix.")
        parser.add_argument("-g",
                            "--gene_info",
                            type=str,
                            required=True,
                            help="The path to the gene info matrix.")
        parser.add_argument("-s",
                            "--samples",
                            type=str,
                            required=False,
                            default=None,
                            help="The path to the expression samples.")
        parser.add_argument("-of",
                            "--outfolder",
                            type=str,
                            required=False,
                            default="output",
                            help="The name of the output folder.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        samples = None
        if self.samples_path is not None:
            print("Loading samples data.")
            samples_df = self.load_file(self.samples_path, header=0, index_col=None)
            samples = set(samples_df.iloc[:, 0].values.tolist())
            del samples_df

        # print("Loading data.")
        # counts_df_list = []
        # stats_df_list = []
        # translate_data = []
        # for filepath in glob.glob(os.path.join(self.data_path, "*")):
        #     df = self.load_file(filepath, header=0, index_col=0)
        #     old_names = list(df.columns)
        #     dataset = None
        #
        #     filename = os.path.basename(filepath)
        #     if "TargetALS" in filename:
        #         df.columns = [re.sub("-", "_", colname) for colname in df.columns]
        #         df.columns = [re.sub("\\.", "_", colname) for colname in df.columns]
        #         df.columns = [re.search(".*(HRA_[0-9]+)", colname).group(1) for colname in df.columns]
        #         dataset = "TargetALS"
        #     elif "Braineac" in filename:
        #         df.columns = [re.search(".*(A653.*)", colname).group(1) for colname in df.columns]
        #         dataset = "Braineac"
        #     elif "GTEx" in filename:
        #         df.columns = [re.search("(.*)_.*", colname).group(1) for colname in df.columns]
        #         dataset = "GTEx"
        #     elif "NABEC" in filename:
        #         df.columns = [re.search(".*_(.*)", colname).group(1) for colname in df.columns]
        #         dataset = "NABEC"
        #     elif "ENA" in filename:
        #         df.columns = [re.search(".*_(.*)", colname).group(1) for colname in df.columns]
        #         dataset = "ENA"
        #     elif "BrainGVEx" in filename:
        #         df.columns = [re.sub("_", "-", colname) for colname in df.columns]
        #         df.columns = [re.sub("\\.", "-", colname) for colname in df.columns]
        #         df.columns = [re.sub("^X", "", colname) for colname in df.columns]
        #         df.columns = [re.search("^([0-9]+-[0-9]+)-[0-9]+-[0-9]+", colname).group(1) for colname in df.columns]
        #         dataset = "BrainGVEx"
        #     elif "BipSeq" in filename:
        #         df.columns = [re.search("Br[0-9]+_(R[0-9]+)", colname).group(1) for colname in df.columns]
        #         dataset = "BipSeq"
        #     elif "UCLA_ASD" in filename:
        #         df.columns = [re.search("[aA-zZ]+[0-9]+_(.+)", colname).group(1) for colname in df.columns]
        #         dataset = "UCLA_ASD"
        #     elif "CMC_HBCC" in filename:
        #         df.columns = [re.search("individualID.*_specimenID.(.*)", colname).group(1) for colname in df.columns]
        #         dataset = "CMC_HBCC"
        #     elif "CMC" in filename:
        #         df.columns = [re.search("^CMC_[aA-zZ]+_[0-9]+_(.*)", colname).group(1) for colname in df.columns]
        #         dataset = "CMC"
        #     elif "MSBB" in filename:
        #         df.columns = [re.sub(".accepted_hits.sort.coordReadsPerGene.out.tab", "", colname) for colname in df.columns]
        #         df.columns = [re.search("AMPAD_MSSM_[0-9]+_(.*)", colname).group(1) for colname in df.columns]
        #         dataset = "MSBB"
        #     elif "ROSMAP" in filename:
        #         df.columns = [re.sub("^X", "", colname) for colname in df.columns]
        #         df.columns = [re.sub("ReadsPerGene.out.tab", "", colname) for colname in df.columns]
        #         df.columns = [re.search(".*_.*_(.*_.*)", colname).group(1) for colname in df.columns]
        #         dataset = "ROSMAP"
        #     elif "MayoCBE" in filename:
        #         df.columns = [re.sub("^X", "", colname) for colname in df.columns]
        #         df.columns = [re.search("[0-9]+_CER", colname).group(0) for colname in df.columns]
        #         dataset = "Mayo_CBE"
        #     elif "MayoTCX" in filename:
        #         df.columns = [re.sub("^X", "", colname) for colname in df.columns]
        #         df.columns = [re.search("[0-9]+_TCX", colname).group(0) for colname in df.columns]
        #         dataset = "Mayo_TCX"
        #     elif "Brainseq" in filename:
        #         # these did not get adjusted in the other talbes, so keep the same
        #         dataset = "Brainseq"
        #     else:
        #         print("Unexpected input file.")
        #         exit()
        #
        #     new_names = list(df.columns)
        #     translate_data.extend(list(zip([dataset] * df.shape[1], old_names, new_names)))
        #
        #     if samples is not None:
        #         print("\tOverlap: {} / {}".format(len(samples.intersection(set(df.columns))), df.shape[1]))
        #
        #     stats_df_list.append(df.iloc[:4, :])
        #     counts_df_list.append(df.iloc[4:, :])
        # counts_df = pd.concat(counts_df_list, axis=1)
        # print(data_df)
        # stats_df = pd.concat(stats_df_list, axis=1).T
        # print(stats_df)
        # translate_df = pd.DataFrame(translate_data, columns=["Dataset", "Original ID", "Revised ID"])
        # print(translate_df)
        #
        # print("Saving data")
        # self.save_file(df=counts_df, outpath=os.path.join(self.file_outdir, "geneCounts.txt.gz"))
        # self.save_file(df=stats_df, outpath=os.path.join(self.file_outdir, "geneCounts_stats.txt.gz"))
        # self.save_file(df=translate_df, outpath=os.path.join(self.file_outdir, "translate.txt.gz"), index=None)
        # del stats_df, translate_df

        print("Loading data")
        counts_df = self.load_file(os.path.join(self.file_outdir, "geneCounts.txt.gz"), header=0, index_col=0)
        translate_df = self.load_file(os.path.join(self.file_outdir, "translate.txt.gz"), header=0, index_col=None)
        std = dict(zip(translate_df.loc[:, "Revised ID"], translate_df.loc[:, "Dataset"]))

        print("Step 1: remove probes with zero variance.")
        mask = counts_df.std(axis=1) != 0
        print("\tUsing {}/{} probes.".format(np.sum(mask), np.size(mask)))
        counts_df = counts_df.loc[mask, :]

        print("Step 2: remove samples with zero counts / variance.")
        mask = (counts_df.std(axis=0) != 0) & (counts_df.sum(axis=0) != 0)
        print("\tUsing {}/{} samples.".format(np.sum(mask), np.size(mask)))
        counts_df = counts_df.loc[:, mask]

        print("Step 3: PCA analysis.")
        self.pca(df=counts_df,
                 sample_to_dataset=std,
                 plot_appendix="_geneCounts")

        print("Step 4: Loading gene length")
        gene_info_df = self.load_file(self.gene_info_path, header=0, index_col=0)
        print(gene_info_df)
        print([gene for gene in counts_df.index if gene not in gene_info_df.index])
        exit()

        print("Step 5: Calculating TPM values")
        # https://btep.ccr.cancer.gov/question/faq/what-is-the-difference-between-rpkm-fpkm-and-tpm/
        # Divide the read counts by the length of each gene in kilobases. This
        # gives you reads per kilobase (RPK).
        kilo_bases_s = gene_info_df.loc[counts_df.index, "MergedExonLength"] / 1e3
        rpk_df = counts_df.divide(kilo_bases_s, axis=0)

        # Count up all the RPK values in a sample and divide this number by
        # 1,000,000. This is your “per million” scaling factor.
        pm_scaling_factor = rpk_df.sum(axis=0) / 1e6

        # Divide the RPK values by the “per million” scaling factor.
        # This gives you TPM.
        tpm_df = rpk_df.divide(pm_scaling_factor, axis=1)

        print("\tSaving data")
        self.save_file(df=tpm_df, outpath=os.path.join(self.file_outdir, "geneCounts.TPM.txt.gz"))

        print("Step 6: PCA analysis.")
        self.pca(df=counts_df,
                 sample_to_dataset=std,
                 plot_appendix="_geneCounts_TPM")

    @staticmethod
    def load_file(inpath, header, index_col, sep="\t", low_memory=True,
                  nrows=None, skiprows=None):
        df = pd.read_csv(inpath, sep=sep, header=header, index_col=index_col,
                         low_memory=low_memory, nrows=nrows, skiprows=skiprows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(inpath),
                                      df.shape))
        return df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, sep="\t"):
        compression = 'infer'
        if outpath.endswith('.gz'):
            compression = 'gzip'

        df.to_csv(outpath, sep=sep, index=index, header=header,
                  compression=compression)
        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def pca(self, df, sample_to_dataset, plot_appendix=""):
        # samples should be on the columns and genes on the rows.
        zscores = (df - df.mean(axis=0)) / df.std(axis=0)
        pca = PCA(n_components=2)
        pca.fit(zscores)
        expl_variance = {"PC{}".format(i+1): pca.explained_variance_ratio_[i] * 100 for i in range(2)}
        components_df = pd.DataFrame(pca.components_)
        components_df.index = ["Comp{}".format(i + 1) for i, _ in enumerate(components_df.index)]
        components_df.columns = df.columns

        print("\tPlotting PCA")
        plot_df = components_df.T
        plot_df["cohort"] = plot_df.index.map(sample_to_dataset)
        self.plot(df=plot_df, x="Comp1", y="Comp2", hue="cohort", palette=self.palette,
                  xlabel="PC1 [{:.2f}%]".format(expl_variance["PC1"]),
                  ylabel="PC2 [{:.2f}%]".format(expl_variance["PC2"]),
                  title="PCA - eigenvectors",
                  filename="eigenvectors_plot{}".format(plot_appendix))

    def plot(self, df, x="x", y="y", hue=None, palette=None, xlabel=None,
             ylabel=None, title="", filename="PCA_plot"):
        if xlabel is None:
            xlabel = x
        if ylabel is None:
            ylabel = y

        sns.set(rc={'figure.figsize': (12, 9)})
        sns.set_style("ticks")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, gridspec_kw={"width_ratios": [0.9, 0.1]})
        sns.despine(fig=fig, ax=ax1)
        ax2.axis('off')

        sns.scatterplot(x=x,
                        y=y,
                        hue=hue,
                        data=df,
                        s=100,
                        linewidth=0,
                        legend=None,
                        palette=palette,
                        ax=ax1)

        ax1.set_title(title,
                      fontsize=20,
                      fontweight='bold')
        ax1.set_ylabel(ylabel,
                       fontsize=14,
                       fontweight='bold')
        ax1.set_xlabel(xlabel,
                       fontsize=14,
                       fontweight='bold')

        if palette is not None:
            handles = []
            for label, color in palette.items():
                if label in df[hue].values.tolist():
                    handles.append(mpatches.Patch(color=color, label=label))
            ax2.legend(handles=handles, loc="center")

        #fig.savefig(os.path.join(self.plot_outdir, "{}.pdf".format(filename)))
        fig.savefig(os.path.join(self.plot_outdir, "{}.png".format(filename)))
        plt.close()

    def print_arguments(self):
        print("Arguments:")
        print("  > Data path: {}".format(self.data_path))
        print("  > Gene info path: {}".format(self.data_path))
        print("  > Samples path: {}".format(self.samples_path))
        print("  > Plot output directory {}".format(self.plot_outdir))
        print("  > File output directory {}".format(self.file_outdir))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
