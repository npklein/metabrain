# TODO: this code still extracts the cohort info from the path instead of from the samplesheet info, has to be
library(data.table)
source('Utility.R')
library(stringr)

combine_QC_metrics <- function(opt){
  AMP_AD_QC <- parse_amp_ad_QC(opt)
  
  combinedMultiQC_data <- parse_multiQC(opt)
  
  STAR_data <- parse_STAR(opt)
  fastQC_data <- parse_fastQC(opt)
  
  combinedMultiQC_data <- get_rnaseq_id_from_SampleFull(combinedMultiQC_data)
  
  # TODO: 1 sample in STAR but not in combinedMultiQC_data: CMC_MSSM_024_MSSM_RNA_BP_PFC_9, figure out why
  multiMetrics_STAR <- merge(combinedMultiQC_data, STAR_data, by=c('SampleFull','cohort'))
  multiMetrics_STAR_AMP_AD <- rbind.fill(multiMetrics_STAR, AMP_AD_QC)
  

  all_pheno_covariate_data <- merge_fastQC_and_phenoData(multiMetrics_STAR_AMP_AD, fastQC_data)
  return(all_pheno_covariate_data)
}

parse_amp_ad_QC <- function(opt){
  AMP_AD_multiMetrics_files <- list.files(path=opt$QcBaseDir,
                                          pattern = 'CombinedMetrics.csv$', 
                                          recursive = T,
                                          full.names = T)
  AMP_AD_multimetrics_qc_all <- data.frame()
  for(f in AMP_AD_multiMetrics_files){
    multimetrics_qc <- fread(f)
    multimetrics_qc$cohort <- get_cohort_from_path(f)
    AMP_AD_multimetrics_qc_all <- rbind(AMP_AD_multimetrics_qc_all, multimetrics_qc)
    AMP_AD_multiMetrics_files <- get_cohort_from_path(f)
  }
  colnames(AMP_AD_multimetrics_qc_all) <- gsub('.*__','',colnames(AMP_AD_multimetrics_qc_all))
  
  colnames(AMP_AD_multimetrics_qc_all)[colnames(AMP_AD_multimetrics_qc_all)=='sample'] <- 'rnaseq_id'
  colnames(AMP_AD_multimetrics_qc_all)[!colnames(AMP_AD_multimetrics_qc_all) %in% c('rnaseq_id','cohort')] <- paste0(colnames(AMP_AD_multimetrics_qc_all)[!colnames(AMP_AD_multimetrics_qc_all) %in% c('rnaseq_id','cohort')],'_alignment')
  return(AMP_AD_multimetrics_qc_all)
}

parse_fastQC <- function(opt){
  print('Readin fastqc files')
  FastQC_multiQC_files <- list.files(path=opt$QcBaseDir, 
                                     pattern = 'multiqc_fastqc.txt$', 
                                     recursive = TRUE,
                                     full.names=T)
  
  FastQC_qc_all <- data.frame()
  for(f in FastQC_multiQC_files){
    FastQC_qc <- fread(f)
    FastQC_qc$cohort <- get_cohort_from_path(f)
    FastQC_qc_all <- rbind(FastQC_qc_all, FastQC_qc, fill=T)
  }
  
  # Because FastQC gives results per fastq file instead of per sample, adjust the names and merge them together (so that paired/single end fastq files are included)
  FastQC_qc_all$Sample <- gsub('.cram','',FastQC_qc_all$Sample)
  
  FastQC_qc_all_R1 <- FastQC_qc_all[grepl('\\.R1$|_1$|_R1$', FastQC_qc_all$Sample),]
  FastQC_qc_all_R1$Sample <- gsub('_1$','', FastQC_qc_all_R1$Sample)
  FastQC_qc_all_R1$Sample <- gsub('\\.R1$','',FastQC_qc_all_R1$Sample)
  FastQC_qc_all_R1$Sample <- gsub('_R1$','',FastQC_qc_all_R1$Sample)
  
  
  FastQC_qc_all_R2 <- FastQC_qc_all[grepl('\\.R2$|*_2$|_R2$', FastQC_qc_all$Sample),]
  FastQC_qc_all_R2$Sample <- gsub('_2$','',FastQC_qc_all_R2$Sample)
  FastQC_qc_all_R2$Sample <- gsub('\\.R2$','',FastQC_qc_all_R2$Sample)
  FastQC_qc_all_R2$Sample <- gsub('_R2$','',FastQC_qc_all_R2$Sample)
  
  # everything that is leftover is single end, save as R1
  single_end <- FastQC_qc_all[!grepl('\\.R1$|_1$|_R1$', FastQC_qc_all$Sample),]
  single_end <- single_end[!grepl('\\.R2$|*_2$|_R2$', single_end$Sample),]
  FastQC_qc_all_R1 <- rbind(FastQC_qc_all_R1, single_end)
  
  
  FastQC_qc_all_merged <- merge(FastQC_qc_all_R1, FastQC_qc_all_R2, by=c('Sample','cohort'), suffix = c('_R1', '_R2'), fill.x=T)
  same_rows(FastQC_qc_all_merged, list(FastQC_qc_all_R1))
  # TargetALS sample naming is different between FastQC and picard results, change them here to make them match later
  if('TargetALS' %in% FastQC_qc_all_merged$cohort){
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample <- str_match( FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample, 
                                                                                         "(CGND-HRA[-|_][0-9]+-2|CGND-HRA[-|_][0-9]+)")[, 2]
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample <-  gsub('-2$','.2',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample)
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample <-  gsub('-','_',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample)
  }
  FastQC_qc_all_merged[,c('File type_R1','File type_R2','Filename_R1','Filename_R2','Encoding_R1','Encoding_R2') := NULL]

  
  return(FastQC_qc_all_merged)
}

parse_STAR <- function(opt){
  print('Reading star QC files')
  STAR_multiQC_files <- list.files(path=opt$QcBaseDir, 
                                   pattern = 'multiqc_star.txt$', 
                                   recursive = TRUE,
                                   full.names = T)
  
  STAR_qc_all <- data.frame()
  for(f in STAR_multiQC_files){
    STAR_qc <- fread(f)
    STAR_qc$cohort <- get_cohort_from_path(f)
    
    STAR_qc_all <- rbind(STAR_qc_all, STAR_qc)
  }
  STAR_qc_all <- data.frame(STAR_qc_all)
  STAR_qc_all$Sample <- gsub('.cram','',STAR_qc_all$Sample)
  colnames(STAR_qc_all)[colnames(STAR_qc_all)=='Sample'] <- 'SampleFull'
  STAR_qc_all[STAR_qc_all$cohort=='CMC',]$SampleFull <- clean_CMC_ID(STAR_qc_all[STAR_qc_all$cohort=='CMC',]$SampleFull)
  return(STAR_qc_all)
}

parse_multiQC <- function(opt){
  RNAseq_metrics <- merge_multiQC_files(opt$QcBaseDir, 'multiqc_picard_RnaSeqMetrics.txt$')
  Multiple_metrics <- merge_multiQC_files(opt$QcBaseDir, 'multiqc_picard_AlignmentSummaryMetrics.txt$')
  print('Merging RNAseqMetrics and MultiMetrics data')
  combined_metrics <- merge(RNAseq_metrics, Multiple_metrics, by=c('SampleFull','cohort',
                                                                   'PF_ALIGNED_BASES_alignment',
                                                                   'PF_ALIGNED_BASES_genotyping'))
  combined_metrics[combined_metrics$cohort=='CMC',]$SampleFull <- clean_CMC_ID(combined_metrics[combined_metrics$cohort=='CMC',]$SampleFull)
  combined_metrics[combined_metrics$cohort=='CMC_HBCC',]$SampleFull <- clean_CMC_ID(combined_metrics[combined_metrics$cohort=='CMC_HBCC',]$SampleFull)
  #same_rows(combined_metrics, list(RNAseq_metrics, Multiple_metrics))
  return(combined_metrics)
}

merge_multiQC_files <- function(path, pattern){
  # Merge multiQC files together
  
  print(paste('merge multiQC_files found in',path, 'with pattern', pattern))
  # All the multiqc files have the same name, so search for the right pattern of files.
  multiQC_files <- list.files(path=path, 
                              pattern = pattern, 
                              recursive = TRUE,
                              full.names = T)
  
  
  print(paste('Parsing',length(multiQC_files),'files'))
  # ENA has 2 metrics: after alignment and after merging BAM files for the same individual before doing genotypeCalling
  # Because of this, we need two version of all the RNAseq metrics columns, one for after alignment, one fore for the genotyping part
  # Want to make sure that _alignment is in the loop first because then we have the fullSample name and can match for the _genotype samples,
  # which lack this information
  multiQC_files_ordered <- c(multiQC_files[grepl('_alignment', multiQC_files)],
                             multiQC_files[grepl('_genotyping', multiQC_files)])
  multiQC_files_ordered <- c(multiQC_files_ordered, multiQC_files[!multiQC_files %in% multiQC_files_ordered])
  
  multiQC_alignment <- data.frame()
  multiQC_genotyping <- data.frame()
  for(f in multiQC_files_ordered){
    cohort <- get_cohort_from_path(f)
    
    multiQC <- fread(f)
    multiQC$cohort <- cohort
    
    # SampleFull is the individual_sample ID
    if(sum('tmp.cram' %in% multiQC$Sample)>0){
      print(f)
    }
    multiQC$SampleFull <- gsub('.cram','',multiQC$Sample)
    
    # Remove columns that are not needed in the output
    multiQC$SAMPLE <- NULL
    multiQC$Sample <- NULL
    multiQC$READ_GROUP <- NULL
    multiQC$LIBRARY <- NULL
    
    # As mentioned before the loop, change the colum names so that there is a difference between QC and genotype QC
    if(cohort == 'ENA'){
      # For ENA the _aligmment and _genotype files need to be merged, but don't know which one is first. So check if sample names are already in the 
      # large table and if so merge the columns together
      if(grepl('multiQC_alignment',  f)) { 
        multiQC <- change_names(multiQC, c('cohort', 'SampleFull'),'_alignment')
      } else if(grepl('multiQC_genotyping',  f)) { 
        multiQC <- change_names(multiQC, c('cohort', 'SampleFull'),'_genotyping')
        multiQC_genotyping <- rbind(multiQC_genotyping, multiQC)        
        # go to next file
        next
      }else{
        stop('should either have multiQC_alignment or multiQC_genotyping in name')
      }
    }else{
      multiQC <- change_names(multiQC, c('cohort', 'SampleFull'),'_alignment')
    }
    
    multiQC_alignment <- rbind.fill(multiQC_alignment, multiQC)
  }
  multiQC_all <- merge(multiQC_alignment, multiQC_genotyping, all.x=T)
  
  same_rows(multiQC_all, list(multiQC_alignment))
  
  # every value for _alignemnt columns should be filled in 
  # (except for PCT_RIBOSOMAL_BASES_alignment , these can be NA in certain situations)
  # double check!
  can_be_na <- c('PCT_RIBOSOMAL_BASES_alignment')
  for(column in names(multiQC_all)[grepl('_alignment',names(multiQC_all))]){
    if(column %in% can_be_na) {next}
    number_of_NA <- sum(is.na(multiQC_all[,column]))
    if(number_of_NA > 0){
      print(paste('Column', column,'is NA for',number_of_NA,'samples'))
      print('head() of those samples:')
      print(head(multiQC_all[which(is.na(multiQC_all[,column])),]))
      t <- multiQC_all[which(is.na(multiQC_all[,column])),]
      stop(paste('ERROR:',column,'should not contain any NA'))
    }
  }
  print('Done!')
  
  return(multiQC_all)
}

merge_fastQC_and_phenoData <- function(multiMetrics_STAR_AMP_AD,fastQC_data){
  ###### merge FastQC and picard metrics together #####

  colnames(fastQC_data)[colnames(fastQC_data)=='Sample'] <- 'rnaseq_id'
  
  all_pheno_covariate_data <- merge(multiMetrics_STAR_AMP_AD, fastQC_data, by=c('rnaseq_id','cohort'), all = TRUE)
  return(all_pheno_covariate_data)
}

get_rnaseq_id_from_SampleFull <- function(combinedMultiQC_data){
  combinedMultiQC_data$rnaseq_id <- NA
  combinedMultiQC_data[combinedMultiQC_data$cohort=="Braineac",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="Braineac",]$Sample, ".*(A653.*)")[, 2]
  combinedMultiQC_data[combinedMultiQC_data$cohort=="GTEx",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="GTEx",]$Sample, "(.*)_.*")[, 2]
  combinedMultiQC_data[combinedMultiQC_data$cohort=="ENA",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="ENA",]$Sample,".*_(.*)")[, 2]
  combinedMultiQC_data[combinedMultiQC_data$cohort=="NABEC",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="NABEC",]$Sample,".*_(.*)")[, 2]
  
  combinedMultiQC_data[combinedMultiQC_data$cohort=="TargetALS",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="TargetALS",]$Sample, ".*(CGND[_|-]HRA.*)")[, 2]
  combinedMultiQC_data[combinedMultiQC_data$cohort=="TargetALS",]$rnaseq_id <- gsub('-2$','.2',combinedMultiQC_data[combinedMultiQC_data$cohort=="TargetALS",]$rnaseq_id)
  combinedMultiQC_data[combinedMultiQC_data$cohort=="TargetALS",]$rnaseq_id <- gsub('-','_',combinedMultiQC_data[combinedMultiQC_data$cohort=="TargetALS",]$rnaseq_id)
  
  
  combinedMultiQC_data[combinedMultiQC_data$cohort=="Brainseq",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="Brainseq",]$Sample, ".*_(R[0-9]+)")[, 2]
  combinedMultiQC_data[combinedMultiQC_data$cohort=="BipSeq",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="BipSeq",]$Sample, ".*_(R[0-9]+)")[, 2]
  combinedMultiQC_data[combinedMultiQC_data$cohort=="CMC",]$rnaseq_id <- str_match(combinedMultiQC_data[combinedMultiQC_data$cohort=="CMC",]$Sample, ".*_([A-Z]+_RNA_.*_[0-9]+$)")[, 2]
  
  CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample <- gsub('_','-',CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample )
  CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample,"^([0-9]+-[0-9]+)-[0-9]+-[0-9]+")[, 2]
  CombinedMetrics[CombinedMetrics$cohort=="UCLA_ASD",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="UCLA_ASD",]$Sample,"[aA-zZ]+[0-9]+_(.+)")[, 2]
  CombinedMetrics[CombinedMetrics$cohort=="CMC_HBCC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="CMC_HBCC",]$Sample,"individualID.*_specimenID.(.*)")[, 2]
    
  return(combinedMultiQC_data)
}
