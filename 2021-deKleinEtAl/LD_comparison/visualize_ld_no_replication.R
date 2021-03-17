#!/usr/bin/Rscript
library(optparse)
library(ggplot2)
library(data.table)
library(reshape2)
library(grid)
library(gtable)
library(ggpubr)
library('scales')


option_list = list(
  make_option(c("-l", "--ld_scores"), type="character",
              help="File containing the LD scores between the top SNPs for each gene of metaBrain that is also found in eqtlGen (fdr < 0.05)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#### Set the data: file to read ####

# This contains the LD scores between the top SNPs for each gene of metaBrain that is also found in eqtlGen (fdr < 0.05)
ld_scores_file <- opt$ld_scores
ld_scores_file <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/LD_comparison/data/metaBrain_eqtlGen_topSNP_LD.txt'

#####

plot_ld <- function(ld_scores){
  # because error messages were in the column they were considered characters, convert to numeric
  ld_scores$`D'` <- as.numeric(ld_scores$`D'`)
  ld_scores$R2 <- as.numeric(ld_scores$R2)
  
  # make the plot
  p <- ggplot(ld_scores, aes(x=`D'`, y=R2, colour=concordant))+
    geom_point(alpha=0.5)+
    xlab("D'")+
    theme_pubr(base_size=18)+
    scale_colour_brewer(palette='Set1')

  
  write.table(ld_scores_merged, file='q',sep='\t','quote'=F, row.names = F)

  ld_scores_merged_disconcordant <- ld_scores_merged[ld_scores_merged$concordant==F & ld_scores_merged$R2 == 1,]
  write.table(ld_scores_merged_disconcordant, file='LD_scores_topSNP_metaBrain_eqtlGen.disconcordant.txt',sep='\t','quote'=F, row.names = F)
  
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  plot_ld(ld_scores_file)
}

