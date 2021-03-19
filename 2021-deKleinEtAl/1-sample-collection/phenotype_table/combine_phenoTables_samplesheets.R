#!/usr/bin/env Rscript

# Get all the different covariate/samplesheet/phenotype files and merge them all together
# while harmonizing column names
options(warn=1)
library(data.table)
library(plyr)
library(optparse)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
pwd <- getwd()
setwd(script.basename)
source('parse_samplesheets.R')
source('utility.R')
source('column_info.R')
setwd(pwd)

test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  opt <- list()
  input_dir <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2020-01-31-alignmentQC-technicalCovariates'
  opt$QcBaseDir <- paste0(input_dir,'/2020-01-31-alignmentQC/')
  opt$metadataDir <- paste0(input_dir,'/2020-02-02-all-cohort-metadata/')
  opt$outputdir <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2020-02-03-phenotype-table/"
  opt$sample_IDs <- paste0('/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2020-01-31-expression-tables/2020-01-31-step4-sample-removal-by-PCA/2020-01-31-step2-second-round-PCA-on-expression/2020-02-10-second-round-includedExpressionSamples.txt')
  return(opt)
}


##### main program functions ######
main <- function(){
  opt <- get_options()
  # (un)comment depending on if it is an interactive test run
  #opt <- test()
  
  all_pheno_covariate_data <- combine_all_data(opt)
  if(length(table(all_pheno_covariate_data$rnaseq_id)[table(all_pheno_covariate_data$rnaseq_id)>1])>0){
    print(table(all_pheno_covariate_data$rnaseq_id)[table(all_pheno_covariate_data$rnaseq_id)>1])
    stop("Some IDs are in double, has to be fixed")
  }
  if(sum(is.na(all_pheno_covariate_data$MetaCohort)>0)){
    stop("MetaCohort is NA for some samples, this should not happen")
  }
  write_output(opt$outputdir, all_pheno_covariate_data)
}


get_options <- function(){
  # Get command line arguments 
  option_list = list(
    make_option(c('-r', '--QcBaseDir'), type='character', default=getwd(), 
                help='Path to RNAseq multiQC data (base directory that contains all subdirs)', 
                metavar='character'),
    make_option(c('-o', '--outputdir'), type='character', default=getwd(), 
                help='Directory to write output to', 
                metavar='character'),
    make_option(c('-s','--sample_IDs'), type='character',
                  help='File with newline separated list of samples that should be included in the phenotype table'),
    make_option(c('-m', '--metadataDir'), type='character', default=getwd(), 
                help='Directory containing metadata of all cohorts', 
                metavar='character')
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt)
}



combine_all_data <- function(opt){
  # first read in samplesheet data. This should theoretically contain all info on IDs that can be used to parse later files
  all_samplesheet_data <- parse_samplesheet_files(opt)
  all_samplesheet_data <- change_rows_if_na(all_samplesheet_data, 'genotype_id',c('individualID','Individual_ID'),no_delete = c('individualID'))
  all_samplesheet_data <- change_rows_if_na(all_samplesheet_data, 'individualID',c('genotype_id'),no_delete = c('genotype_id'))
  
  check_if_all_samples_included(opt, all_samplesheet_data)
  
  return(all_samplesheet_data)
}
check_if_all_samples_included <- function(opt, all_pheno_covariate_data){
  expression_samples <- fread(opt$sample_IDs,header=F)
  expression_samples$V1 <- gsub('_mergedWithResequenced','',expression_samples$V1)
  s_not_found <- expression_samples[!expression_samples$V1 %in% all_pheno_covariate_data$rnaseq_id]$V1
  

  all_pheno_covariate_data[(grepl('481_120514',all_pheno_covariate_data$rnaseq_id)),]$cohort
  if(length(s_not_found)>0){
    cat(paste0(length(s_not_found)," samples not found"))
    cat("heaD() of samples not found:\n")
    cat(head(s_not_found))
    cat("\n")
    stop("Not all samples were found in QC/covariate table")
  }
  else{cat("All samples found")}
}

write_output <- function(outputdir, combined_metrics){
  # make sure that spaces are removed from column nnames
  colnames(combined_metrics) <- gsub(' ','_', colnames(combined_metrics))
  colnames(combined_metrics) <- gsub('%','perc', colnames(combined_metrics))
  # Rename some columns
  
  combined_metrics <-   Filter(function(x)!all(is.na(x)), combined_metrics)

  
  
  # Change the order of columns so that identifiers are first
  change_order <- c('SampleFull','rnaseq_id','MetaCohort','cohort','BroadBrainRegion','SpecificBrainRegion')
  for(c in change_order){
    setcolorder(combined_metrics, c(c, setdiff(names(combined_metrics), c)))
    
  }
  
  outfile <- paste0(outputdir,'/',Sys.Date(),'.brain.phenotypes.txt')
  write.table(combined_metrics, outfile, quote=F, sep='\t', row.names=F)
  print(paste('Written to',outfile))
  
  column_info <- get_column_info(combined_metrics)
  outfile <- paste0(outputdir,'/',Sys.Date(),'.brain.phenotypes_columnDescription.txt')
  write.table(column_info, outfile, quote=F, sep='\t', row.names=F)
  print(paste('Written to',outfile))
}
#####


##### file parsing functions #####




# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main()
}

