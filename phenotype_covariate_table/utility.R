# Utility functions to parse and merge phenotype and covariate data


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
                metavar='character'),
    make_option(c('-s', '--sampleList'), type='character', default=getwd(), 
                help='File with list of IDs for checking that they are all included in the phenotable', 
                metavar='character')
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  return(opt)
}
test <- function(){
  # This function is only here for development so that I can more easily test interactively in Rstudio
  opt <- list()
  input_dir <- '/Users/NPK/UMCG/projects/biogen/cohorts/'
  opt$QcBaseDir <- input_dir
  opt$metadataDir <- paste0(input_dir,'/joined_analysis/phenotype_QC_and_covariates_table/input/all_cohort_metadata')
  opt$outputdir <- paste0(input_dir,'/joined_analysis/phenotype_QC_and_covariates_table/output/')
  opt$sample_IDs <- paste0(input_dir,'/joined_analysis/phenotype_QC_and_covariates_table/input/2019-04-11-Freeze2.TMM.SampleFilter.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors-pc1to4-nonoutliers-forcovariates.txt')
  return(opt)
}
same_rows <- function(df, list_of_df){
  # helper function testing if all df haave same number of rows
  for(df2 in list_of_df){
    if(nrow(df)!=nrow(df2)){
      stop('ERROR: df should have same number of rows as other df, but did not')
    }
  }
}
get_cohort_from_path <- function(path){
  # order is important, goes from most to leas specific. 
  # e.g. CMC_HBCC needs to be before CMC, otherwise would return CMC due to grep
  if(grepl('ENA',path)) { return('ENA') }
  if(grepl('NABEC',path)) { return('NABEC') }
  if(grepl('GTEx',path)) { return('GTEx') }
  if(grepl('TargetALS',path)) { return('TargetALS') }
  if(grepl('Braineac',path)) { return('Braineac') }
  if(grepl('Brainseq',path)) { return('Brainseq') }
  if(grepl('CMC_HBCC',path)) { return('CMC_HBCC') }
  if(grepl('CMC',path)) { return('CMC') }
  if(grepl('BipSeq',path)) { return('BipSeq') }
  if(grepl('BrainGVEx',path)) { return('BrainGVEx') }
  if(grepl('UCLA_ASD',path)) { return('UCLA_ASD') }
  if(grepl('Mayo_CBE',path)) { return('Mayo_CBE') }
  if(grepl('Mayo_TCX',path)) { return('Mayo_TCX') }
  if(grepl('ROSMAP',path)) { return('ROSMAP') }
  if(grepl('MSSM',path)) { return('MSSM') }
  stop(paste('unknown cohort:',path))
}
clean_CMC_ID <- function(CMC_IDs){
  CMC_IDs <- gsub('specimenID.','',gsub('individualID.','',CMC_IDs))
  return(CMC_IDs)
}
na_count <- function(df){
  count <- sapply(df, function(y) sum(length(which(is.na(y)))))
  return(count)
}
change_names <- function(dt, c, a){
  # quick helper function for changing all names except for a few listed in a vector
  # dt: data table to change names of
  # c: vector with column names which NOT to change
  # a: what has to be appended to the columns
  names(dt)[which(!names(dt) %in% c)] <- paste0(names(dt)[which(!names(dt) %in% c)], a)
  return(dt)
}
parse_rnaseq_ID_from_fullName <- function(){
  # tmp function, hopefully do not need it and can instead parse it out of the samplesheet files
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
}

change_rows_if_na <- function(dt, column_to_keep, column_to_merge_into, overwrite=F, no_delete = c()){
  # Check if column_to_keep is NA, if so loop over columns in column_to_merge_into and try to itteratively fill it
  if(!is.data.table(dt)){
    stop('Error: input to change_rows_if_na has to be data.table')
  }
  # first check if not multiple columns are not NA 
  if(nrow(dt[length(column_to_merge_into)-rowSums(is.na(dt[,..column_to_merge_into]))>1,])){
    # multiple columns not NA, have to implement this first
    if(!overwrite){
      stop('multiple column not NA, implement this first')
    }
  }
  
  for(column_name in column_to_merge_into){
    # Get the rows where column_to_keep == NA
    rows_na <- which(is.na(dt[,..column_to_keep]) | dt[,..column_to_keep] == '')
    # Then fill in with the value of the next column in the loop (might be NA, but doesn't matter). Do this for all columns to merge
    dt[rows_na, (column_to_keep) := dt[rows_na,][[column_name]]]
  }
  # Then delete all columns that were merged into the column_to_keep that are not in no_delete
  column_to_delete <- column_to_merge_into[!column_to_merge_into %in% no_delete]
  if(length(column_to_delete) > 0){
    dt[,(column_to_delete):=NULL]
  }
  return(dt)
}

fix_broadman_area <- function(dt){
  dt[BrodmannArea=='Amygdala',BrodmannArea := 'ba25' ]
  dt[, BrodmannArea := gsub('bm','ba',tolower(BrodmannArea))]
  
  dt[organism_part == 'Brodmann (1909) area 24',BrodmannArea := 'ba24']
  dt[organism_part == 'Brodmann (1909) area 9', BrodmannArea := 'ba9']
  dt[body_site_annotations == 'Brain - Anterior cingulate cortex (BA24)',BrodmannArea := 'ba24']
  dt[body_site_annotations == 'Brain - Frontal Cortex (BA9)', BrodmannArea := 'ba9']
  dt[organism_part == 'Brodmann (1909) area 24', organism_part := NA ]
  dt[organism_part == 'Brodmann (1909) area 9', organism_part := NA ]
  dt[body_site_annotations == 'Brain - Anterior cingulate cortex (BA24)',body_site_annotations := NA]
  dt[body_site_annotations == 'Brain - Frontal Cortex (BA9)', body_site_annotations := NA]
  t <- apply(dt, 1, function(r) any(r %in% c("Brain_-_Frontal_Cortex_(BA9)")))
  return(dt)
}
  
convert_to_broad_region <- function(dt){
  dt[,SpecificBrainRegion := gsub(' ','_',SpecificBrainRegion)]
  
  dt[SpecificBrainRegion %in% c('DLPFC','dorsolateral_prefrontal_cortex'), SpecificBrainRegion := 'Dorsolateral_Prefrontal_Cortex']
  dt[SpecificBrainRegion %in% c('Cerebral_Frontal_Cortex'), SpecificBrainRegion := 'Prefrontal_Cortex']
  dt[SpecificBrainRegion %in% c('temporal_cortex','TemporalCortex'), SpecificBrainRegion := 'Temporal_Cortex']
  dt[SpecificBrainRegion %in% c('amygdala'), SpecificBrainRegion := 'Amygdala']
  dt[SpecificBrainRegion %in% c('cerebellum'), SpecificBrainRegion := 'Cerebellum']
  dt[SpecificBrainRegion %in% c('Cortex_Frontal','frontal_cortex'), SpecificBrainRegion := 'Cerebral_Frontal_Cortex']
  dt[SpecificBrainRegion %in% c('HIPP'), SpecificBrainRegion := 'hippocampus']
  dt[SpecificBrainRegion %in% c('PUTM'), SpecificBrainRegion := 'putamen']
  dt[SpecificBrainRegion %in% c('SNIG'), SpecificBrainRegion := 'substantia_nigra']
  dt[SpecificBrainRegion %in% c('sACC'), SpecificBrainRegion := 'subgenual_anterior_cingulate_cortex']
  
  CX <- c('ba10','ba36','ba44','ba22','ba24', 'ba9',"cerebral_cortex","Cortex","Cortex_Occipital","Cortex_Sensory","Dorsolateral_Prefrontal_Cortex","Motor_Cortex",
          "Motor_Cortex_Lateral","Motor_Cortex_Medial","parietal_cortex","Prefrontal_Cortex","Temporal_Cortex","visual_cortex",'subgenual_anterior_cingulate_cortex')
  
  CER <- c('Cerebellum','cerebellar_hemisphere')
  LIMBIC_SYSTEM <- c('hippocampus','hippocampus_proper','hypothalamus','Amygdala')
  BASAL_GANGLIA <- c('substantia_nigra','caudate_nucleus','putamen','PUTM_SNIG','nucleus_accumbens')
  SPINAL_CORD <- c("C1_segment_of_cervical_spinal_cord","Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Spinal_Cord_Thoracic" )
  dt$BroadBrainRegion <- ''
  dt[SpecificBrainRegion %in% CX, BroadBrainRegion := 'CX']
  dt[SpecificBrainRegion %in% CER, BroadBrainRegion := 'CER']
  dt[SpecificBrainRegion %in% LIMBIC_SYSTEM, BroadBrainRegion := 'LIMBIC_SYSTEM']
  dt[SpecificBrainRegion %in% BASAL_GANGLIA, BroadBrainRegion := 'BASAL_GANGLIA']
  dt[SpecificBrainRegion %in% SPINAL_CORD, BroadBrainRegion := 'SPINAL_CORD']
  dt[SpecificBrainRegion=='Liver', BroadBrainRegion:='Liver']
  dt[grepl('Cell_Line', SpecificBrainRegion), BroadBrainRegion:='Cell line']
  dt[BroadBrainRegion=='',BroadBrainRegion:=NA]
  
  
  return(dt)
}
