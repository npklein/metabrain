

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
  options <- test()
  PCs <- annotate_PCs(options)
  
  
  plot_PCAs(PCs, options)
}

get_options <- function(){
  # Get command line arguments 
  option_list = list(
    make_option(c("-p", "--PCA"), type="character", 
                help="path to file with PCAs (it is assumed to be a table with first column sample name, rest columns PC scores)", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir (needs to exist and contain a figures/ subdir)", metavar="character"),
    make_option(c("-c", "--covariates"), type="character", default=NULL,
                help="Table with technical covariates", 
                metavar="character"),
    make_option(c("-m", "--maxCovCor"), type="character", default=NULL,
                help="Table with maximum correlation values between pca and covariates", 
                metavar="character")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  options = parse_args(opt_parser);
  return(options)
}

test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  options <- list()
  options$PCA <- "pc1_2.txt"
  options$output <- getwd()
  options$covariates <- '2020-02-05-freeze2dot1.TMM.Covariates.withBrainRegion.txt'
  options$maxCovCor <- '2020-02-05-CovariateCorrelations.txt-max.txt'
  return(options)
}

annotate_PCs <- function(options){
  PCs <- fread(options$PCA)
  covariates <- fread(options$covariates)
  
  PCs_with_covariates <- merge(PCs, covariates, by='Sample', all.x=T)
  return(PCs_with_covariates)
}

plot_max_cov(max_cov){
  max_cov <- max_cov[order(-max_cov$MaxRSq),]
  max_cov$included <- F
  max_cov[1:10,]$included <- T
  ggplot(max_cov, aes(reorder(Covariate,MaxRSq),MaxRSq, fill="red", colour='black'))+
    geom_bar(stat='identity', aes(alpha=included))+
    coord_flip()+
    geom_errorbar(aes(ymin=MeanRSq-VarRsq, ymax=MeanRSq+VarRsq),
                  width=.7,                    # Width of the error bars
                  position=position_dodge(.9))+
    scale_fill_manual(values=c("#fee8c8"))+
    scale_colour_manual(values=c("black"))+
    theme_bw(base_size=12)+
    guides(alpha=F, fill=F, colour=FALSE)+
    xlab('')+
    scale_alpha_manual(values=c(0.05,1))
  ggsave('figures/cov_correlation.pdf',height=12, width=6)
}

plot_PCAs <- function(PCs, options){
  
  PCs <- PCs[,colSums(is.na(PCs))<nrow(PCs), with=FALSE]
  PCs <- PCs[!is.na(PCs$cohort),]
  # removing things that are not of interest for correlation (e.g. SampleFull), or that have most samples
  # in one group so correlation wouldn't give any inf (e.g. per_sequence_quality_scores_R2: pass=6970, fail=38, warn=5)
  covariates_only <- PCs[,!c('Sample','Filename_R1', 'Filename_R2','FILTER','pcaFiltered',
                            'File.type_R2','File.type_R1','FastQC_original_sample','SampleFull',
                            'Encoding_R1','Encoding_R2','lib_selection','per_sequence_quality_scores_R2',
                            'per_base_sequence_content_R2','basic_statistics_R2','kmer_content_R2',
                            'per_sequence_quality_scores_R1','basic_statistics_R1','kmer_content_R1',
                            'study')]

  max_cov <- read.table( options$maxCovCor, header=T, sep='\t')
  plot_max_cov(max_cov)
  
  max_cov <- max_cov[order(-max_cov$MaxRSq),]
  max_cov$included <- F
  max_cov[1:10,]$included <- T
  
  plot_list<- list()  
  i <-  1
  for(cov in c('cohort','BroadBrainRegion',as.character(max_cov_dot2[max_cov$included==T,]$Covariate))){
    print(cov)
    covariates_only_sub <- covariates_only[,c('Comp1','Comp2',cov),with=F]
  
    p <- ggplot(covariates_only_sub, aes_string('Comp1', 'Comp2', fill=cov))+
      geom_point(pch=21)+
      theme_bw()
    if(cov != 'cohort' & cov != 'BroadBrainRegion'){
      p <- p+scale_fill_gradientn(colours = rainbow(5))
    }
    plot_list[[i]] <- p
    i <- i+1
  }
  
  library(gridExtra)
  pdf('figures/top10_cov_plus.pdf',width=16,height=12)
  do.call(grid.arrange,c(plot_list, ncol=3))
  dev.off()
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main()
}
