

library(data.table)
library(ggplot2)
library(GGally)
library(ggrepel)
library(optparse)


main <- function(){
  opt <- get_options()
  # (un)comment depending on if it is an interactive test run
  opt <- test()
  PCs <- annotated_PCs(opt)
  plot_PCAs(opt)
}

get_options <- function(){
  # Get command line arguments 
  option_list = list(
    make_option(c("-p", "--PCA"), type="character", default=getwd(), 
                help="path to file with PCAs", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir (needs to exist and contain a figures/ subdir)", metavar="character"),
    make_option(c("-a", "--annotation"), type="character", default=NULL,
                help="Directory that contains annotation files of the samples (e.g. brain region/batch/etc)", 
                metavar="character")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt)
}

test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  opt <- list()
  input_dir <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/expression_PCA/'
  opt$PCA <- paste0(input_dir,"TargetALs_Braineac_CMC_MayoCBE_MayoTCX_MSSM_ROSMAP_GTEx.ENA.TMM.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.PCAOverSamplesEigenvectors.txt.gz")
  opt$output <- input_dir
  opt$annotation <- paste0(input_dir,'2019-03-06-ENA-annotation/')
  return(opt)
}

annotate_PCs <- function(opt){
  PCs <- fread(opt$PCA)

  rna_annotation <- fread(paste0(opt$annotation,'/rna-tissue.annot.txt'), header=F)
  rna_dataset <- fread(paste0(opt$annotation,'/rna-dataset.annot.txt'),header=F)
  #genotype <- fread(paste0(opt$annotation,'rna-genotype.annot.txt'),header=F)
  colnames(rna_annotation) <- c('sample','region')

  PCs <- data.frame(PCs)
  rownames(PCs) <- PCs$X.
  PCs$X. <- NULL

  PCs <- PCs[1:10]
  PCs$region <- factor(rna_annotation[match(rownames(PCs), rna_annotation$sample),]$region)
  PCs$cohort <- factor(rna_dataset[match(rownames(PCs), rna_dataset$V1),]$V2)
  #PCs$genotype <- factor(genotype[match(rownames(PCs), genotype$V1),]$V2)
  return(PCs)
}

make_pca_plot <- function(PC_x, PC_y, name, colour, opt){
  p <- ggplot(PCs, aes_string(x=gsub('PC','Comp',PC_x), 
                              y=gsub('PC','Comp',PC_y), 
                              colour=colour))+
    geom_point()+
    theme_pubr(base_size=20)+
    xlab(PC_x)+
    ylab(PC_y)+
    scale_colour_brewer(palette="Dark2")+ 
    guides(colour = guide_legend(override.aes = list(size=5)))


  pdf(paste0(opt$output,'/figures/',name,'PCA.pdf'), width=12, height=8)
  p <- ggExtra::ggMarginal(p, type = "histogram",  groupFill = TRUE)
  dev.off()
  return(p)
}

plot_PCAs <- function(opt){
  p1 <- make_pca_plot('PC1','PC2','PC1_PC2_tissue','region',opt)
  p2 <- make_pca_plot('PC3','PC4', 'PC3_PC4_tissue','region',opt)
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main()
}