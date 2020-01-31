# Get all the different covariate/samplesheet/phenotype files and merge them all together
# while harmonizing column names
options(warn=2)
library(data.table)
library(plyr)
source('parse_samplesheets.R')
source('utility.R')
source('parse_QC_metrics.R')
source('column_info.R')

##### main program functions ######
main <- function(){
  #opt <- get_options()
  # (un)comment depending on if it is an interactive test run
  opt <- test()
  
  all_pheno_covariate_data <- combine_all_data(opt)
  
  write_output(opt$outputdir, all_pheno_covariate_data)
}

combine_all_data <- function(opt){
  # first read in samplesheet data. This should theoretically contain all info on IDs that can be used to parse later files
  all_samplesheet_data<- parse_samplesheet_files(opt)
  # QC metrics are already in different script, comment out for now and check if is necesarry to add in here later
  #all_qc_metrics <- combine_QC_metrics(opt)
  # remove cohort because this is in the samplesheet data
  #all_qc_metrics$cohort <- NULL
  
  #all_samplesheet_data <- data.frame(all_samplesheet_data)

  #all_pheno_covariate_data <- merge(all_samplesheet_data,all_qc_metrics, by='rnaseq_id',fill.x=T)
  all_pheno_covariate_data <- all_samplesheet_data

  all_pheno_covariate_data <- change_rows_if_na(all_pheno_covariate_data, 'genotype_id',c('individualID','Individual_ID'),no_delete = c('individualID'))
  all_pheno_covariate_data <- change_rows_if_na(all_pheno_covariate_data, 'individualID',c('genotype_id'),no_delete = c('genotype_id'))
  
  all_pheno_covariate_data <- check_and_select_all_samples_included(opt, all_pheno_covariate_data)
  
  return(all_pheno_covariate_data)
}

check_and_select_all_samples_included <- function(opt, dt){
  samples <- fread(opt$sample_IDs,header=F)
  dt <- dt[!duplicated(dt$rnaseq_id),]
  dt <- dt[rnaseq_id %in% samples$V1,]
  if(length(samples[!samples$V1 %in% dt$rnaseq_id]$V1)){
    print('Mssing:')
    print(samples[!samples$V1 %in% dt$rnaseq_id]$V1)
    stop('For some samples phenotype data is missing')
  }

  dt <- fix_genotype_IDs(dt,opt)
  
  return(dt)
}

fix_genotype_IDs <- function(dt,opt){
  dt[individualID=='?', individualID := NA,]
  dt[genotype_id=='?', genotype_id := NA,]
  genotype_files <- list.files(opt$genotype_IDs_dir, pattern='*txt',full.names = T)
  genotype_IDs <- data.frame()
  for(f in genotype_files){
    if(grepl('Integrative-Individuals.txt',f)){
      next
    }
    genotype_IDs_current_cohort <- fread(f, header=F)
    cohort <- gsub('-Individuals.txt','',basename(f))
    genotype_IDs_current_cohort$cohort <- cohort
    type <- 'array'
    if(grepl('WES',f)){
      type='WES'
    }
    genotype_IDs_current_cohort$type <- type
    genotype_IDs <- rbind(genotype_IDs,genotype_IDs_current_cohort )
  }
  genotype_IDs[,V1 := as.character(V1)]
  genotype_IDs <- genotype_IDs[!duplicated(genotype_IDs$V1),]
  integrative_genotype <- fread(paste0(opt$genotype_IDs_dir,'/Integrative-Individuals.txt'),header=F)
  
  # R function
  search = function(x) {
    genotype_IDs_cohort_subset <- genotype_IDs[grepl(unique(x[['MetaCohort']]), genotype_IDs$cohort),]
    if(unique(x[['MetaCohort']])=='Braineac'){
      x[['individualID']] <- paste0(x[['individualID']], '_1')
    }
    
    if(unique(x[['MetaCohort']])=='TargetALS'){
      x[['individualID']] <- paste0(x[['individualID']], '-b38')
    }
    
    pattern <- paste0(x[['individualID']],'$|',gsub('-','_',x[['individualID']]),'$')
    full_id <- grep(pattern, genotype_IDs_cohort_subset[type=='array',]$V1, value=T)
    if(length(full_id) > 1){
      print(full_id)
      print(pattern)
      stop('can not return more than 1')
    }
    
    if(length(full_id)==0){
      
      # try also in WES
      full_id <- grep(pattern, 
                      genotype_IDs_cohort_subset[type=='WES',]$V1, value=T)
      if(length(full_id)==0){
        return(NA)
      }
      if(length(full_id) > 1){
        print(full_id)
        print(pattern)
        stop('can not return more than 1')
      }
    }
    return(full_id)
  }
  dt[MetaCohort=='NABEC' ,genotype_id := unlist(apply(dt[MetaCohort=='NABEC',], 1, search))]
  dt[MetaCohort=='CMC' ,genotype_id := unlist(apply(dt[MetaCohort=='CMC',], 1, search))]
  dt[MetaCohort=='Braineac' ,genotype_id := unlist(apply(dt[MetaCohort=='Braineac',], 1, search))]
  dt[MetaCohort=='TargetALS' & !genotype_id %in%genotype_IDs$V1 ,genotype_id := unlist(apply(dt[MetaCohort=='TargetALS' & !genotype_id %in%genotype_IDs$V1,], 1, search))]
  
  table(dt[!genotype_id %in% genotype_IDs$V1 & ! genotype_id %in% integrative_genotype$V1,]$MetaCohort)
  dt[!genotype_id %in% genotype_IDs$V1,]$genotype_id
  dt[cohort=='BipSeq' & !genotype_id %in% genotype_IDs$V1,]$individualID
  
  
  # add some missing sample IDs that hj found
  missing_samples <- fread(paste0(opt$metadataDir,'/missing_sample_links.txt'))
  
  dt[rnaseq_id == missing_samples$V2, genotype_id := missing_samples$V1]
  dt[rnaseq_id=='Br1023_R3489',]$genotype_id
  
  return(dt)
}

write_output <- function(outputdir, phenotype_data){
  # make sure that spaces are removed from column nnames
  colnames(phenotype_data) <- gsub(' ','_', colnames(phenotype_data))
  colnames(phenotype_data) <- gsub('%','perc', colnames(phenotype_data))
  # Rename some columns
  phenotype_data <- phenotype_data[,which(unlist(lapply(phenotype_data, function(x)!all(is.na(x))))),with=F]

  
  # Change the order of columns so that identifiers are first
  change_order <- c('BrodmannArea','BroadBrainRegion','SpecificBrainRegion','cohort','MetaCohort','SampleFull','genotype_id','rnaseq_id','individualID')
  for(c in change_order){
    setcolorder(phenotype_data, c(c, setdiff(names(phenotype_data), c)))
  }
  
  outfile <- paste0(outputdir,'/2019-07-29-brain.phenotypes.txt')
  write.table(phenotype_data, outfile, quote=F, sep='\t', row.names=F)
  print(paste('Written to',outfile))
  
  column_info <- get_column_info(phenotype_data)
  outfile <- paste0(outputdir,'/brain.phenotypes.columnDescription.txt')
  write.table(column_info, outfile, quote=F, sep='\t', row.names=F)
  print(paste('Written to',outfile))
}
#####


##### file parsing functions #####




# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main(eqtlGen_replication, ld_scores_file, fdr_threshold)
}
