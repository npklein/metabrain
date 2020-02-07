

library(data.table)
library(ggplot2)
library(GGally)
library(ggrepel)
library(optparse)
library(ggExtra)
library(RColorBrewer)

main <- function(){
  options <- get_options()
  # (un)comment depending on if it is an interactive test run
  options <- test()
  PCs <- annotate_PCs(options)
  
  PCs[which(PCs$age_death=='90+'),]$age_death <- '90'
  PCs$age_death <- as.numeric(PCs$age_death)
  PCs$PMI_.in_hours. <- as.numeric(PCs$PMI_.in_hours.)
  PCs$RIN <- as.numeric(PCs$RIN)
  plot_PCAs(PCs, options)
}

get_options <- function(){
  # Get command line arguments 
  option_list = list(
    make_option(c("-p", "--PCA"), type="character", default=getwd(), 
                help="path to file with PCAs", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir (needs to exist and contain a figures/ subdir)", metavar="character"),
    make_option(c("-c", "--covariates"), type="character", default=NULL,
                help="Table with technical covariates", 
                metavar="character"),
    make_option(c("-p", "--phenotype"), type="character", default=NULL,
                help="Phenotype table", 
                metavar="character"),
    make_option(c("-n", "--numberOfPcaToPlot"), type="character", default=2,
                help="The number of PCAs to plot (default: 2)", 
                metavar="integer")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  options = parse_args(opt_parser);
  return(options)
}

test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  options <- list()
  input_dir <- '/Users/NPK/UMCG/git_projects/brain_eQTL/expression_PCA/'
  options$PCA <- paste0(input_dir,"pc1_2.txt")
  options$output <- input_dir
  options$annotation <- paste0(input_dir,'2019-03-06-ENA-annotation/')
  options$phenotype <- '2020-02-03.brain.phenotypes.txt'
  options$numberOfPcaToPlot <- 2
  return(options)
}

annotate_PCs <- function(options){
  PCs <- fread(options$PCA)

  #rna_annotation <- fread(paste0(options$annotation,'/rna-tissue.annot.txt'), header=F)
  #rna_dataset <- fread(paste0(options$annotation,'/rna-dataset.annot.txt'),header=F)
  
  phenotable <- fread(options$phenotype)
  
  #genotype <- fread(paste0(options$annotation,'rna-genotype.annot.txt'),header=F)
  #colnames(rna_annotation) <- c('sample','region')

  PCs <- data.frame(PCs)
  rownames(PCs) <- PCs$X.
  PCs$X. <- NULL

  PCs <- PCs[1:(options$numberOfPcaToPlot+1)]
  #PCs$region <- rna_annotation[match(rownames(PCs), rna_annotation$sample),]$region
  #PCs$cohort <- rna_dataset[match(rownames(PCs), rna_dataset$V1),]$V2
  #PCs[grepl('SRR',rownames(PCs)),]$region <- 'Brain'

  #PCs$region <- factor(PCs$region)
  #PCs$cohort <- factor(PCs$cohort)
  PCs$SampleFull <- rownames(PCs)
  
  PCs_with_pheno <- merge(PCs, phenotable, by='SampleFull', all.x=T)
  #PCs$genotype <- factor(genotype[match(rownames(PCs), genotype$V1),]$V2)
  return(PCs_with_pheno)
}


make_pca_plot <- function(PCs, PC_x, PC_y, column){
  if(is.numeric(PCs[[column]])){
    # add the super small nuber in case y-axis needs to be log scaled, but only for those columns where the max value > 100000 (so that columns with e.g. percentages don't get log scaled)
    if(max(PCs[!is.na(PCs[column]),][column])> 10000){
      PCs[!is.na(PCs[column]),][column] <-PCs[!is.na(PCs[column]),][column] +0.00000000000000000000001
    }
    PCs[[column]] <- as.numeric(PCs[[column]])
  }else{
    colourCount = length(unique(PCs[[column]]))
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
    dark2ExtendedPallete <- getPalette(colourCount)
  }
  
  p <- ggplot(PCs, aes_string(x=gsub('PC','Comp',PC_x), 
                         y=gsub('PC','Comp',PC_y), 
                         colour=column))+
    geom_point(alpha=0.5)+
    theme_bw(base_size=20)+
    xlab(PC_x)+
    ylab(PC_y)
  if(is.numeric(PCs[[column]])){
    p <- p+ scale_color_distiller(palette = "RdPu",na.value="grey80")
    if(max(PCs[!is.na(PCs[column]),][column])> 10000){
      p <- p +  guides(colour=guide_legend(paste0('log(',column,')')))
    }
  }else{
    p <- p+scale_colour_manual(values=dark2ExtendedPallete,na.value="grey80")+
      guides(colour = guide_legend(override.aes = list(size=5)))
  }

  

  return(p)
}

plot_PCAs <- function(PCs, options){
  
  PCs <- PCs[,colSums(is.na(PCs))<nrow(PCs)]
  
  PCs_subset <- PCs[12:length(colnames(PCs))]
  PCs_subset <- PCs_subset[!colnames(PCs_subset) %in% c('individualIdentifier.inferred','genotype_id','projid','Flowcell','study','msex',
                                                        'cohort','rnaseq_id')]
  columns_in_front <- c('MEDIAN_5PRIME_BIAS_alignment', 'MEDIAN_3PRIME_BIAS_alignment', 'PCT_PF_READS_ALIGNED_alignment','total_reads')
  for(c in columns_in_front){
    if(c %in% colnames(PCs_subset)){
      PCs_subset <- PCs_subset[c(c, setdiff(names(PCs_subset), c))]
    }
    
  }
  
  pdf(paste0(options$output,'figures/all_pca.pdf'),width=12, height=8, onefile = TRUE)
  # read-length, read-depth, percentage mapped reads 
  for(column in colnames(PCs_subset)){
    if(nrow(PCs)/sum(is.na(PCs[[column]])) < 1){
      print(paste('skip',column))
      next
    }
    if(column == ''){
      next
    }
    name <- paste0('PC1_PC_2_',column)
    outfile <- paste0(options$output,'/figures/all_pca/',name,'.PCA.png')
    print(paste(column,'-',outfile))
    
    
    p <- make_pca_plot(PCs, 'PC1', 'PC2', column)
    print(p)

    ggsave(outfile, plot=p, width=12, height=8)
    #p <- ggMarginal(p, type = "histogram",  groupFill = TRUE)
    #pdf(outfile, width=12, height=8)
    #print(p)
    #dev.off()
    #print(paste('written to',outfile))
  }
  dev.off()
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main()
}