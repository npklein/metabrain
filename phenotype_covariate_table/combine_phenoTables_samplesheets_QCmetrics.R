# Get all the different covariate/samplesheet/phenotype files and merge them all together
# while harmonizing column names
options(warn=2)
library(data.table)
library(plyr)

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
  #opt <- get_options()
  # (un)comment depending on if it is an interactive test run
  opt <- test()
  
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
  combinedMultiQC_data <- parse_multiQC(opt)
  STAR_data <- parse_STAR(opt)
  fastQC_data <- parse_fastQC(opt)
  
  #TODO: fix
  amp_ad <- parse_amp_ad(opt)
  amp_ad$SampleFull <- amp_ad$rnaseq_id
  
  
  all_pheno_covariate_data <- merge(combinedMultiQC_data, STAR_data, by=c('SampleFull','cohort'))
  all_pheno_covariate_data <- rbind.fill(all_pheno_covariate_data, amp_ad)
  
  all_pheno_covariate_data$rnaseq_id <- NA
  #TODO: fix the covariates data
  #all_covariates <- parse_covariate_files(opt)
  #head(STAR_data[which(!STAR_data$SampleFull %in% all_covariates$SampleFull),c('SampleFull','cohort')])
  
  all_pheno_covariate_data <- merge_fastQC_and_phenoData(all_pheno_covariate_data, fastQC_data)
  
  
  #all_pheno_covariate_data <- merge(all_pheno_covariate_data, all_covariates, by=c('SampleFull','cohort'))
  #same_rows(all_pheno_covariate_data, list(combinedMultiQC_data, STAR_data))
  return(all_pheno_covariate_data)
}

write_output <- function(outputdir, combined_metrics){
  # make sure that spaces are removed from column nnames
  colnames(combined_metrics) <- gsub(' ','_', colnames(combined_metrics))
  colnames(combined_metrics) <- gsub('%','perc', colnames(combined_metrics))
  # Rename some columns
  colnames(combined_metrics)[colnames(combined_metrics)=='TotalReads'] <- 'total_reads'
  colnames(combined_metrics)[colnames(combined_metrics)=='TOTAL_READS_alignment'] <- 'total_reads'
  
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
merge_fastQC_and_phenoData <- function(all_pheno_covariate_data,fastQC_data){
  ###### merge FastQC and picard metrics together #####
  #all_pheno_covariate_data[grepl('PCT', colnames(all_pheno_covariate_data))] <- all_pheno_covariate_data[grepl('PCT', 
  #                                                                                     colnames(all_pheno_covariate_data))]/100

  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='Braineac',]$rnaseq_id <- str_match(all_pheno_covariate_data[all_pheno_covariate_data$cohort=='Braineac',]$SampleFull, '.*(A653.*)')[, 2]
  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='GTEx',]$rnaseq_id <- str_match(all_pheno_covariate_data[all_pheno_covariate_data$cohort=='GTEx',]$SampleFull, '(.*)_.*')[, 2]
  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='ENA',]$rnaseq_id <- str_match(
                all_pheno_covariate_data[all_pheno_covariate_data$cohort=='ENA',]$SampleFull,'.*_(.*)')[, 2]
  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='NABEC',]$rnaseq_id <- str_match(
                all_pheno_covariate_data[all_pheno_covariate_data$cohort=='NABEC',]$SampleFull,'.*_(.*)')[, 2]
  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='Brainseq',]$rnaseq_id <- str_match(
                all_pheno_covariate_data[all_pheno_covariate_data$cohort=='Brainseq',]$SampleFull,'.*_(.*)')[, 2]
  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='TargetALS',]$rnaseq_id <- str_match(all_pheno_covariate_data[all_pheno_covariate_data$cohort=='TargetALS',]$SampleFull, '.*(HRA[_|-][0-9]+)')[, 2]
  all_pheno_covariate_data[all_pheno_covariate_data$cohort=='CMC',]$rnaseq_id <- str_match(all_pheno_covariate_data[all_pheno_covariate_data$cohort=='CMC',]$SampleFull, '.*_([a-zA-Z]+_RNA.*)')[, 2]

  colnames(FastQC_qc_all_merged)[colnames(FastQC_qc_all_merged)=='Sample'] <- 'rnaseq_id'
  FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$rnaseq_id <- str_match(FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$rnaseq_id, '.*(HRA[_|-][0-9]+)')[, 2]
  
    
  all_pheno_covariate_data <- merge(all_pheno_covariate_data, FastQC_qc_all_merged, by=c('rnaseq_id','cohort'), all = TRUE)
  return(all_pheno_covariate_data)
}

parse_covariate_files <- function(opt){
  # AMP_AD
  AMP_AD <- parse_amp_ad(opt)
  
  Braineac <- fread(paste0(opt$metadataDir,'/Braineac/SampleInfoBiogen.csv'))
  colnames(Braineac)[colnames(Braineac)=='Sample_ID'] <- 'rnaseq_id'
  colnames(Braineac)[colnames(Braineac)=='SD No'] <- 'genotype_id'
  Braineac$genotype_id <- gsub('/','_', Braineac$genotype_id)
  Braineac$cohort <- 'Braineac'
  
  Brainseq <- fread(paste0(opt$metadataDir,'/Brainseq/phenotypeFile_LIBD_szControl.csv'))
  colnames(Brainseq)[colnames(Brainseq)=='RNum'] <- 'rnaseq_id'
  colnames(Brainseq)[colnames(Brainseq)=='BrNum'] <- 'genotype_id'
  Brainseq$cohort <- 'Brainseq'
  
  CMC <- parse_CMC(opt)
  CMC$cohort <- 'CMC'
  
  GTEx <- fread(paste0(opt$metadataDir,'/GTEx/GTEx-RNASeq-samplesheet.txt'))
  colnames(GTEx)[colnames(GTEx)=='Comment[ENA_RUN]'] <- 'rnaseq_id'
  colnames(GTEx)[colnames(GTEx)=='Characteristics[individual]'] <- 'genotype_id'
  GTEx$cohort <- 'GTEx'
  
  NABEC <- parse_NABEC(opt)
  NABEC$cohort <- 'NABEC'
  
  TargetALS <- parse_TargetALS(opt)
  TargetALS$cohort <- 'TargetALS'
  
  #TODO: fix below
  merged_covariates <- data.frame(AMP_AD)
  for(df in list(data.frame(Brainseq), data.frame(CMC), data.frame(GTEx), data.frame(NABEC), data.frame(TargetALS),data.frame(Braineac))){
    # Because rbind.fill was giving trouble I now only rbind columns present in all datasets, fix this later
    merged_covariates <- merged_covariates[which(colnames(merged_covariates) %in% colnames(df))]
    df <- df[which(colnames(df) %in% colnames(merged_covariates))]
    merged_covariates <- rbind.fill(merged_covariates, df)
  }
  
  merged_covariates$SampleFull <- paste0(merged_covariates$genotype_id,'_',merged_covariates$rnaseq_id)
  
  return(merged_covariates)
}

parse_ROSMAP <- function(opt){
  ROSMAP_IDs <- fread(paste0(opt$metadataDir,'/AMP_AD/ROSMAP_IDkey.csv'), na.strings=c(""))
  # Only need thos with an RNAseq ID
  ROSMAP_IDs <- ROSMAP_IDs[!is.na(ROSMAP_IDs$rnaseq_id),]
  # select informative columns
  ROSMAP_IDs <- ROSMAP_IDs[,c('projid','wgs_id','rnaseq_id')]
  
  ROSMAP_clinical <- fread(paste0(opt$metadataDir,'/AMP_AD/ROSMAP_clinical.csv'))
  ROSMAP <- merge(ROSMAP_IDs, ROSMAP_clinical, by='projid', all.x=T)
  same_rows(ROSMAP, list(ROSMAP_IDs))
  colnames(ROSMAP)[colnames(ROSMAP)=='wgs_id'] <- 'genotype_id'

  return(ROSMAP)
}

parse_Mayo <- function(opt){
  Mayo_wgs <- fread(paste0(opt$metadataDir,'/AMP_AD/AMP-AD_Mayo_WGS_covariates.csv'))
  Mayo_RNAseq_CBE <- fread(paste0(opt$metadataDir,'/AMP_AD/MayoRNAseq_RNAseq_CBE_covariates.csv'))
  Mayo_RNAseq_TCX <- fread(paste0(opt$metadataDir,'/AMP_AD/MayoRNAseq_RNAseq_TCX_covariates.csv'))
  colnames(Mayo_RNAseq_TCX)[colnames(Mayo_RNAseq_TCX)=='ID'] <- 'rnaseq_id'
  colnames(Mayo_RNAseq_CBE)[colnames(Mayo_RNAseq_CBE)=='SampleID'] <- 'rnaseq_id'
  colnames(Mayo_RNAseq_TCX)[colnames(Mayo_RNAseq_TCX)=='Gender'] <- 'Sex'
  colnames(Mayo_RNAseq_TCX)[colnames(Mayo_RNAseq_TCX)=='FLOWCELL'] <- 'Flowcell'
  if(sum(Mayo_RNAseq_CBE$SampleID %in% Mayo_RNAseq_TCX$SampleID) > 0){ stop ('Should only be different SampleID in CBE and TCX files')}
  Mayo_RNAseq <- rbind(Mayo_RNAseq_CBE, Mayo_RNAseq_TCX)
  
  # Have to merge the wgs with the RNAseq data, change column names and remove doubles (some are in both wgs and RNAseq table)
  colnames(Mayo_wgs)[colnames(Mayo_wgs)=='RNAseq_CER_SampleID'] <- 'rnaseq_id'
  Mayo_CER <- merge(Mayo_wgs, Mayo_RNAseq_CBE, by=c('rnaseq_id', 'Diagnosis','AgeAtDeath','ApoE','PMI','Sex'))
  
  Mayo_wgs$rnaseq_id <- NULL
  colnames(Mayo_wgs)[colnames(Mayo_wgs)=='RNAseq_TCX_SampleID'] <- 'rnaseq_id'
  Mayo_TCX<- merge(Mayo_wgs, Mayo_RNAseq_TCX, by=c('rnaseq_id', 'Diagnosis','AgeAtDeath','ApoE','PMI','Sex'))
  
  Mayo_CER$RNAseq_TCX_SampleID <- NULL
  Mayo <- rbind(Mayo_CER, Mayo_TCX)
  for(c in colnames(Mayo)[grepl('TCX',colnames(Mayo))]){
    Mayo[[c]] <- NULL
    Mayo[[gsub('TCX','CER', c)]] <- NULL
  }
  
  colnames(Mayo)[colnames(Mayo)=='WGS_Participant_ID'] <- 'genotype_id'
  
  return(Mayo)
}

parse_MSBB <- function(opt){
  MSBBB_covariates <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_WES_covariatestsv'))
  colnames(MSBBB_covariates)[colnames(MSBBB_covariates)=='Region'] <- 'BrodmannArea'
  MSBB_update <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_RNAseq_covariates_November2018Update.csv'))
  MSBB_update <- MSBB_update[MSBB_update$fileType=='fastq',]
  MSBB_update$synapseId <- NULL
  MSBB_update$fileName <- NULL
  MSBB_update$fileType <- NULL
  if(nrow(MSBB_update) != length(unique(MSBB_update$sampleIdentifier))) { stop("should only be unique samples")}
  
  MSBB_covariates <- rbind.fill(MSBB_update, MSBBB_covariates)
  MSBB_covariates$BrodmannArea <- gsub('-','',MSBB$BrodmannArea)
  
  MSBB_clinical <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_clinical.csv'))
  
  MSBB <- merge(MSBB_covariates, MSBB_clinical, by='individualIdentifier',all.x=T)
  
  same_rows(MSBB, list(MSBB_covariates))
  colnames(MSBB)[colnames(MSBB)=='individualIdentifier'] <- 'genotype_id'
  colnames(MSBB)[colnames(MSBB)=='sampleIdentifier'] <- 'rnaseq_id'
  return(MSBB)
}

parse_amp_ad_QC <- function(opt){
  AMP_AD_multiMetrics_files <- list.files(path=opt$QcBaseDir,
                                          pattern = 'CombinedMetrics.csv$', 
                                          recursive = T,
                                          full.names = T)
  AMP_AD_multimetrics_qc_all <- data.frame()
  for(f in AMP_AD_multiMetrics_files){
    multimetrics_qc <- fread(f)
    multimetrics_qc$cohort <- sapply(strsplit(basename(f), '_CombinedMetrics'), '[[', 1)
    AMP_AD_multimetrics_qc_all <- rbind(AMP_AD_multimetrics_qc_all, multimetrics_qc)
  }
  colnames(AMP_AD_multimetrics_qc_all) <- gsub('.*__','',colnames(AMP_AD_multimetrics_qc_all))
  names(AMP_AD_multimetrics_qc_all)[names(AMP_AD_multimetrics_qc_all) == 'sample'] <- 'Sample'
    
  colnames(AMP_AD_multimetrics_qc_all)[colnames(AMP_AD_multimetrics_qc_all)=='Sample'] <- 'rnaseq_id'
  colnames(AMP_AD_multimetrics_qc_all)[!colnames(AMP_AD_multimetrics_qc_all) %in% c('rnaseq_id','cohort')] <- paste0(colnames(AMP_AD_multimetrics_qc_all)[!colnames(AMP_AD_multimetrics_qc_all) %in% c('rnaseq_id','cohort')],'_alignment')
  return(AMP_AD_multimetrics_qc_all)
}

parse_amp_ad <- function(opt){
  # TODO: harmonize all columns
  Mayo <- parse_Mayo(opt)
  Mayo$cohort <- 'Mayo'
  ROSMAP <- parse_ROSMAP(opt)
  ROSMAP$cohort <- 'ROSMAP'
  MSBB <- parse_MSBB(opt)
  MSBB$cohort <- 'MSBB'
  AMP_AD <- rbind.fill(Mayo, ROSMAP)
  AMP_AD <- rbind.fill(AMP_AD, MSBB)
  AMP_AD_QC <- parse_amp_ad_QC(opt)
  
  AMP_AD <- merge(AMP_AD, AMP_AD_QC, by=c('rnaseq_id','cohort'),
                  all=T)
  
  return(AMP_AD)
}

parse_CMC <- function(opt){
  CMC_clinical <- fread(paste0(opt$metadataDir,'/CMC/CMC-CMC_HBCC_clinical_.csv'))
  colnames(CMC_clinical) <- gsub(' ','_',colnames(CMC_clinical))
  CMC_rnaseq <- fread(paste0(opt$metadataDir,'/CMC/CMC_Human_DLPFC_rnaSeq.csv'))
  CMC_snp <- fread(paste0(opt$metadataDir,'/CMC/CMC_Human_SNP.csv'))
  CMC_dna_metadata <- fread(paste0(opt$metadataDir,'/CMC/CMC_MSSM-Penn-Pitt_DLPFC_DNA-metaData.csv'))
  CMC_rna_metadata <- fread(paste0(opt$metadataDir,'/CMC/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv'))
  
  colnames(CMC_rna_metadata)[colnames(CMC_rna_metadata)=='DLPFC_RNA_Sequencing_Sample_ID'] <- 'rnaseq_id'
  colnames(CMC_rnaseq)[colnames(CMC_rnaseq)=='Sample_RNA_ID'] <- 'rnaseq_id'
  
  
  CMC_rnaseq <- rbind(CMC_rna_metadata[!CMC_rna_metadata$rnaseq_id %in% CMC_rnaseq$rnaseq_id,])
  
  if(nrow(CMC_snp)!=length(unique(CMC_snp$Individual_ID)) | nrow(CMC_snp)!=length(unique(CMC_snp$Genotyping_Sample_ID))){ stop("More than 1 DNA per individual, change code")}
  CMC <- merge(CMC_rnaseq, CMC_snp, by=c('Individual_ID'), all.x=T)
  
  CMC <- merge(CMC, CMC_clinical, by='Individual_ID',all.x=T)
  CMC <- merge(CMC, CMC_dna_metadata, by='Individual_ID',all.x=T)
  colnames(CMC)[colnames(CMC)=='Individual_ID'] <- 'genotype_id'
  same_rows(CMC, list(CMC_rnaseq))
  return(CMC)
}

parse_NABEC <- function(opt){
  NABEC <- fread(paste0(opt$metadataDir,'/NABEC/NABEC_phenotypes.txt'))
  

  NABEC_RNA <- NABEC[NABEC$Assay_Type=='RNA-Seq',]
  NABEC_RNA$wgs_id <- gsub('-[a-z]+-fctx.*','',NABEC_RNA$Sample_Name)
  NABEC_RNA$rnaseq_id <- NABEC_RNA$wgs_id
  NABEC_WXS <- NABEC[NABEC$Assay_Type=='WXS',]
  NABEC_WXS$wgs_id <- gsub('_e','',gsub('-wes.*','',NABEC_WXS$Sample_Name))
  
  NABEC <- merge(NABEC_RNA, NABEC_WXS, by='wgs_id')
  colnames(NABEC)[colnames(NABEC)=='wgs_id'] <- 'genotype_id'
  return(NABEC)
}

parse_TargetALS <- function(opt){
  TargetALS_metadata_rna <- fread(paste0(opt$metadataDir,'/TargetALS/RNA_Metadata_TALS_2_5July2018.txt'))
  TargetALS_metadata_wgs <- fread(paste0(opt$metadataDir,'/TargetALS/TALS_Data_Release_WGS_Metadata_10Oct2018.txt'))
  TargetALS_clinical1 <- fread(paste0(opt$metadataDir,'/TargetALS/UMG-Target_ALS_RNA_Clinical_Data_06122018.txt'))
  TargetALS_clinical2 <- fread(paste0(opt$metadataDir,'/TargetALS/UMG-Target_ALS_WGS_Clinical_Data_06122018.txt'))
  
  # TODO: merge on all overlapping columns
  TargetALS <- merge(TargetALS_metadata_rna,TargetALS_metadata_wgs, by='ExternalSubjectId',all.x=T)
  #TODO: add clinical data
  #TargetALS <- merge(TargetALS,TargetALS_clinical1, by=colnames(TargetALS)[which(colnames(TargetALS) %in% colnames(TargetALS_clinical1))],all.x=T)
  #TargetALS <- merge(TargetALS,TargetALS_clinical2, by=colnames(TargetALS)[which(colnames(TargetALS) %in% colnames(TargetALS_clinical2))],all.x=T)
  
  colnames(TargetALS)[colnames(TargetALS)=='ExternalSampleId'] <- 'rnaseq_id'
  colnames(TargetALS)[colnames(TargetALS)=='ExternalSubjectId'] <- 'genotype_id'
  return(TargetALS)
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
  
  FastQC_qc_all_R1 <- FastQC_qc_all[grepl('*\\.R1$|_1$|_R1$', FastQC_qc_all$Sample),]
  FastQC_qc_all_R1$Sample <- gsub('_1$','', FastQC_qc_all_R1$Sample)
  FastQC_qc_all_R1$Sample <- gsub('\\.R1$','',FastQC_qc_all_R1$Sample)
  FastQC_qc_all_R1$Sample <- gsub('_R1$','',FastQC_qc_all_R1$Sample)
  
  
  FastQC_qc_all_R2 <- FastQC_qc_all[grepl('*\\.R2$|*_2$|_R2$', FastQC_qc_all$Sample),]
  FastQC_qc_all_R2$Sample <- gsub('_2$','',FastQC_qc_all_R2$Sample)
  FastQC_qc_all_R2$Sample <- gsub('\\.R2$','',FastQC_qc_all_R2$Sample)
  FastQC_qc_all_R2$Sample <- gsub('_R2$','',FastQC_qc_all_R2$Sample)
  
  
  FastQC_qc_all_merged <- merge(FastQC_qc_all_R1, FastQC_qc_all_R2, by=c('Sample','cohort'), suffix = c('_R1', '_R2'), fill.x=T)
  same_rows(FastQC_qc_all_merged, list(FastQC_qc_all_R1))
  # TargetALS sample naming is different between FastQC and picard results, change them here to make them match later
  if('TargetALS' %in% FastQC_qc_all_merged$cohort){
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample <- str_match(FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=='TargetALS',]$Sample , '.*(HRA[_|-][0-9]+)')[, 2]
  }
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
    
    change_names <- function(dt, c, a){
      # quick helper function for changing all names except for a few listed in a vector
      # dt: data table to change names of
      # c: vector with column names which NOT to change
      # a: what has to be appended to the columns
      names(dt)[which(!names(dt) %in% c)] <- paste0(names(dt)[which(!names(dt) %in% c)], a)
      return(dt)
    }
    
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

get_column_info <- function(combined_metrics){
  # Add info on what each of the columns mean, and check that there is info included for all columns in the output file
  add_info <- function(df, col1, col2){
    df <- rbind(df, data.frame(col1,col2))
    return(df)
  }
  info <- data.frame(col1=c(), col2=c())
  info <- add_info(info, 'SampleFull','Sample name as used in the Compute5 pipeline we ran. Often in form of <individual ID>_<Sample/run_ID>, but not always')
  info <- add_info(info, 'cohort','Name of the cohort. This is the larger cohort, e.g. cohort PsychEncode consists of studies BipSeq, BrainVEX and others, and cohort ENA consists of many studies')
  for (c in c('PCT_UTR_BASES','UTR_BASES','PCT_CODING_BASES','RIBOSOMAL_BASES','PF_NOT_ALIGNED_BASES','PCT_MRNA_BASES','CORRECT_STRAND_READS','MEDIAN_CV_COVERAGE','PCT_RIBOSOMAL_BASES','CODING_BASES','PCT_INTRONIC_BASES','MEDIAN_3PRIME_BIAS','PF_ALIGNED_BASES','NUM_R1_TRANSCRIPT_STRAND_READS','MEDIAN_5PRIME_BIAS','PCT_USABLE_BASES','PF_BASES','NUM_UNEXPLAINED_READS','IGNORED_READS','PCT_R2_TRANSCRIPT_STRAND_READS','PCT_R1_TRANSCRIPT_STRAND_READS','INTERGENIC_BASES','INTRONIC_BASES','NUM_R2_TRANSCRIPT_STRAND_READS','MEDIAN_5PRIME_TO_3PRIME_BIAS','PCT_INTERGENIC_BASES','PCT_CORRECT_STRAND_READS','INCORRECT_STRAND_READS')){
  
    info <- add_info(info, paste0(c,'_alignment'), 'QC after alignment - From PicardTools RNAseqMetrics, see https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics for details')
    info <- add_info(info, paste0(c,'_genotyping'), 'QC after MarkDuplicates - From PicardTools RNAseqMetrics, see https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics for details')
    
  }
  for(c in c('PF_HQ_ALIGNED_READS','READS_ALIGNED_IN_PAIRS','MEAN_READ_LENGTH','BAD_CYCLES','PF_INDEL_RATE','PCT_READS_ALIGNED_IN_PAIRS','STRAND_BALANCE','PF_HQ_ALIGNED_Q20_BASES','CATEGORY','TOTAL_READS','PCT_CHIMERAS','PF_HQ_ERROR_RATE','PF_ALIGNED_BASES','PCT_PF_READS','PF_READS_IMPROPER_PAIRS','PF_READS','PCT_ADAPTER','PCT_PF_READS_ALIGNED','PF_NOISE_READS','PF_HQ_MEDIAN_MISMATCHES','PF_HQ_ALIGNED_BASES','PCT_PF_READS_IMPROPER_PAIRS','PF_MISMATCH_RATE','PF_READS_ALIGNED')){
    info <- add_info(info, paste0(c,'_alignment'), 'QC after alignment - From different PicardMetrics, see https://broadinstitute.github.io/picard/picard-metric-definitions.html for details')
    info <- add_info(info, paste0(c,'_genotyping'), 'QC after MarkDuplicates - From different PicardMetrics, see https://broadinstitute.github.io/picard/picard-metric-definitions.html for details')
  }
  
  colnames(info) <- c('Column name', 'Description')
  
  for(column in colnames(combined_metrics)){
    if(! column %in% info$`Column name`){
      stop(paste('No description of column',column))
    }
  }
  # Make order in info same as in results
  info <- info[match(colnames(combined_metrics), info$`Column name`),]
  return(info)
}

not_included_yet <- function(opt){


  
  all_cohorts <- all_cohorts[order(all_cohorts$cohort),]
  all_cohorts$sampleID <- 1:nrow(all_cohorts)
  #######
  
  ###### set filter based on 3 measurements ######
  all_cohorts$FILTER <- 'NO'
  print('set filter')
  if(sum(!is.na(all_cohorts$PCT_CODING_BASES) & all_cohorts$PCT_CODING_BASES <= 0.1) > 0){
    all_cohorts[!is.na(all_cohorts$PCT_CODING_BASES) & all_cohorts$PCT_CODING_BASES <= 0.1,]$FILTER <- 'YES'
  }
  if(sum(!is.na(all_cohorts$PCT_PF_READS_ALIGNED) & all_cohorts$PCT_PF_READS_ALIGNED <= 0.6) > 0){
    all_cohorts[!is.na(all_cohorts$PCT_PF_READS_ALIGNED) & all_cohorts$PCT_PF_READS_ALIGNED <= 0.6,]$FILTER <- 'YES'
  }
  if(sum(!is.na(all_cohorts$uniquely_mapped_percent) & all_cohorts$uniquely_mapped_percent <= 60) > 0){
    all_cohorts[!is.na(all_cohorts$uniquely_mapped_percent) & all_cohorts$uniquely_mapped_percent <= 60,]$FILTER <- 'YES'
  }
  print('done')
  
  all_cohorts[is.na(all_cohorts$SampleFull),]$SampleFull <- all_cohorts[is.na(all_cohorts$SampleFull),]$Sample
  merged_metrics_files <- paste0(opt$output,'/all_cohort_STAR_RNAseqMetrics_MultipleMetrics_FastQC.txt')

  nSamples <- all_cohorts %>% group_by(cohort) %>% dplyr::summarize(no = sum(FILTER=='NO'),yes = sum(FILTER=='YES'))
  nSamples$total <- nSamples$no + nSamples$yes
  
  #####
}
#####

##### Utility functions #####
same_rows <- function(df, list_of_df){
  # helper function testing if all df haave same number of rows
  for(df2 in list_of_df){
    if(nrow(df)!=nrow(df2)){
      stop('ERROR: df should have same number of rows as other df, but did not')
    }
  }
}

get_cohort_from_path <- function(path){
  if(grepl('ENA',path)) { return('ENA') }
  if(grepl('NABEC',path)) { return('NABEC') }
  if(grepl('GTEx',path)) { return('GTEx') }
  if(grepl('TargetALS',path)) { return('TargetALS') }
  if(grepl('Braineac',path)) { return('Braineac') }
  if(grepl('Brainseq',path)) { return('Brainseq') }
  if(grepl('CMC',path)) { return('CMC') }
  stop(paste('unknown cohort:',path))
}

clean_CMC_ID <- function(CMC_IDs){
  CMC_IDs <- gsub('specimenID.','',gsub('individualID.','',CMC_IDs))
  return(CMC_IDs)
}
#####

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main(eqtlGen_replication, ld_scores_file, fdr_threshold)
}
