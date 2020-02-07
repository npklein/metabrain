library(data.table)
library(ggplot2)
library(GGally)
library(ggrepel)
library(optparse)
library(ggExtra)
library(RColorBrewer)
library(fastDummies)
library("psych")

main <- function(){
  options <- get_options()
  # (un)comment depending on if it is an interactive test run
#  options <- test()
  PCs <- annotate_PCs(options)
  
  plot_PCAs(PCs, options)
}

get_options <- function(){
  # Get command line arguments 
  option_list = list(
    make_option(c("-p", "--PCA"), type="character", 
                help="path to file with PCAs (it is assumed to be a table with first column sample name, rest columns PC scores)"), 
    make_option(c("-c", "--covariates"), type="character", 
                help="Table with technical covariates"),
    make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir (needs to exist and contain a figures/ subdir)")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  options = parse_args(opt_parser);
  return(options)
}

test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  options <- list()
  options$PCA <- "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-sample-removal-by-PCA/2020-01-31-step3-third-round-PCA-on-expression/pc1_2.txt"
  options$covariates <- "/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-01-31-alignmentQC-technicalCovariates/2020-02-04-freeze2dot1.TMM.Covariates.txt"
  return(options)
}

annotate_PCs <- function(options){
  PCs <- fread(options$PCA)
  covariates <- fread(options$covariates)
  
  PCs_with_covariates <- merge(PCs, covariates, by='Sample', all.x=T)
  return(PCs_with_covariates)
}


plot_PCAs <- function(PCs, options){
  
  PCs <- PCs[,colSums(is.na(PCs))<nrow(PCs), with=FALSE]
  PCs <- PCs[!is.na(PCs$cohort),]
  # removing things that are not of interest for correlation (e.g. SampleFull), or that have most samples
  # in one group so correlation wouldn't give any inf (e.g. per_sequence_quality_scores_R2: pass=6970, fail=38, warn=5)
  covariates_only <- PCs[,!c('Sample','Comp1','Comp2','Filename_R1', 'Filename_R2','FILTER','pcaFiltered',
                            'File.type_R2','File.type_R1','FastQC_original_sample','SampleFull',
                            'Encoding_R1','Encoding_R2','lib_selection','per_sequence_quality_scores_R2',
                            'per_base_sequence_content_R2','basic_statistics_R2','kmer_content_R2',
                            'per_sequence_quality_scores_R1','basic_statistics_R1','kmer_content_R1')]
  
 
  covariates_only_with_dummies <- dummy_cols(covariates_only)
  numeric <- which(sapply(covariates_only_with_dummies,is.numeric))
  covariates_only_with_dummies_numeric <- covariates_only_with_dummies[,numeric,with=F]
  covariates_only_with_dummies_numeric <- covariates_only_with_dummies_numeric[!is.na(covariates_only_with_dummies_numeric$cohort_GTEx),]
  
  PC1_correlated <- corr.test(PCs$Comp1, covariates_only_with_dummies_numeric,use='complete.obs')
  PC2_correlated <- corr.test(PCs$Comp2, covariates_only_with_dummies_numeric,use='complete.obs')


  PC1_r_melt <- melt(PC1_correlated$r)
  PC1_p_melt <- melt(PC1_correlated$p)
  PC1 <-  merge(PC1_r_melt, PC1_p_melt, by='Var2')
  PC1$Var1.x <- NULL
  PC1$Var1.y <- NULL
  PC1 <- PC1[!(is.na(PC1$value.x) | is.na(PC1$value.y)),]
  PC1$PC <- 1

  PC2_r_melt <- melt(PC2_correlated$r)
  PC2_p_melt <- melt(PC2_correlated$p)
  PC2 <-  merge(PC2_r_melt, PC2_p_melt, by='Var2')
  PC2$Var1.x <- NULL
  PC2$Var1.y <- NULL
  PC2 <- PC2[!(is.na(PC2$value.x) | is.na(PC2$value.y)),]
  PC2$PC <- 2
  print(head(PC1))
  print(head(PC2))
  colnames(PC1) <- c('Covariate','R', 'pval','PC')
  colnames(PC2) <- c('Covariate','R', 'pval','PC')
  
  
  cov_plot <- function(pc,pc_num){
    ggplot(pc, aes(x=reorder(Covariate,R),y=R,
                   fill=pval < 0.05,colour='black'))+
      geom_bar(stat='identity', position='dodge')+
      theme_bw(12)+ 
      coord_flip()+
      scale_fill_manual(values=c('#d8b365','#5ab4ac'))+
      scale_colour_manual(values=c('black'))+
      guides(colour = FALSE)+
      xlab('')+
      ylab(paste0('Pearson correlation with PC',pc_num))+ 
      theme(
        legend.position = c(0.2, 0.98),
      )
  }
  p_PC1 <- cov_plot(PC1,1)
  ggsave('figures/correlated_covariates_PC1.pdf',width=6, height=24)
  p_PC2 <- cov_plot(PC2,2)
  ggsave('figures/correlated_covariates_PC2.pdf',width=6, height=24)
    
  
  
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main()
}
