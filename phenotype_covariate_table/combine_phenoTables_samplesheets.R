#!/usr/bin/env Rscript

# Get all the different covariate/samplesheet/phenotype files and merge them all together
# while harmonizing column names
options(warn=2)
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
source('parse_QC_metrics.R')
source('column_info.R')
setwd(pwd)

test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  opt <- list()
  input_dir <- '/Users/NPK/UMCG/projects/biogen/cohorts/'
  opt$QcBaseDir <- input_dir
  opt$metadataDir <- paste0(input_dir,'/joined_analysis/phenotype_QC_and_covariates_table/all_cohort_metadata')
  opt$outputdir <- paste0(input_dir,'/joined_analysis/phenotype_QC_and_covariates_table')
  return(opt)
}


##### main program functions ######
main <- function(){
  opt <- get_options()
  # (un)comment depending on if it is an interactive test run
  #opt <- test()
  
  all_pheno_covariate_data <- combine_all_data(opt)

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
  all_samplesheet_data<- parse_samplesheet_files(opt)
  sum(names(all_samplesheet_data) == 'rnaseq_id')
  # QC metrics are already in different script, comment out for now and check if is necesarry to add in here later
  #all_qc_metrics <- combine_QC_metrics(opt)
  # remove cohort because this is in the samplesheet data
  #all_qc_metrics$cohort <- NULL
  
  #all_samplesheet_data <- data.frame(all_samplesheet_data)

  #all_pheno_covariate_data <- merge(all_samplesheet_data,all_qc_metrics, by='rnaseq_id',fill.x=T)
  all_pheno_covariate_data <- all_samplesheet_data

  change_rows_if_na(all_pheno_covariate_data, 'genotype_id',c('individualID','Individual_ID'),no_delete = c('individualID'))
  change_rows_if_na(all_pheno_covariate_data, 'individualID',c('genotype_id'),no_delete = c('genotype_id'))
  
  check_if_all_samples_included(opt, all_pheno_covariate_data)
  
  return(all_pheno_covariate_data)
}
check_if_all_samples_included <- function(opt, dt){
  samples <- fread(opt$sample_IDs,header=F)
  samples[!samples$V1 %in% all_pheno_covariate_data$rnaseq_id]$V1
  
  all_pheno_covariate_data[grepl('AN04479_BA7',rnaseq_id),]$rnaseq_id
}

write_output <- function(outputdir, combined_metrics){
  # make sure that spaces are removed from column nnames
  colnames(combined_metrics) <- gsub(' ','_', colnames(combined_metrics))
  colnames(combined_metrics) <- gsub('%','perc', colnames(combined_metrics))
  # Rename some columns
  
  combined_metrics <- combined_metrics[, colSums(is.na(combined_metrics)) != nrow(combined_metrics)]
  
  
  # Change the order of columns so that identifiers are first
  change_order <- c('cohort','SampleFull')
  for(c in change_order){
    setcolorder(combined_metrics, c(c, setdiff(names(combined_metrics), c)))
    
  }
  
  outfile <- paste0(outputdir,'/brain.phenotype_QC_covariates.txt')
  write.table(combined_metrics, outfile, quote=F, sep='\t', row.names=F)
  print(paste('Written to',outfile))
  
  column_info <- get_column_info(combined_metrics)
  outfile <- paste0(outputdir,'/brain.phenotype_QC_covariates.columnDescription.txt')
  write.table(column_info, outfile, quote=F, sep='\t', row.names=F)
  print(paste('Written to',outfile))
}
#####


##### file parsing functions #####




# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main()
}

