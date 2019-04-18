# add description of columns in dataframe
get_column_info <- function(combined_metrics){
  number_NA <- na_count(combined_metrics)
  # Add info on what each of the columns mean, and check that there is info included for all columns in the output file
  add_info <- function(df, colname, description, coltype, AMP_AD, Braineac, Brainseq, CMC, ENA, GTEx, NABEC, TargetALS, encoding=''){
    # for each MetaCohort a boolean if that measurement was taken for that MetaCohort
    df <- rbind(df, data.frame(colname,description, coltype, encoding,  number_NA[colname], AMP_AD, Braineac, Brainseq, CMC, ENA, GTEx, NABEC, TargetALS))
      return(df)
  }
  info <- data.frame(col1=c(), col2=c())
  
  # Add sample info column descriptions
  info <- add_info(info, 'rnaseq_id','RNAseq identifier that matches FastQC ID','SampleInfo', 1,1,1,1,1,1,1,1)
  info <- add_info(info, 'SampleFull',
                   'Sample name as used in the Compute5 pipeline we ran. Often in form of <individual ID>_<Sample/run_ID>, but not always',
                   'SampleInfo', 1,1,1,1,1,1,1,1)
  info <- add_info(info, 'cohort','Name of the cohort. This is the larger cohort, e.g. cohort PsychEncode consists of smaller cohorts BipSeq, BrainVEX and others, and cohort ENA consists of many cohorts',
                   'SampleInfo', 1,1,1,1,1,1,1,1)
  info <- add_info(info, 'MetaCohort','Name of the larger cohort, e.g. cohort PsychEncode, which consists out of multiple smaller cohorts. For smaller cohorts, see cohort column',
                   'SampleInfo', 1,1,1,1,1,1,1,1)
  info <- add_info(info, 'Gender','Reported gender of the sample. For genotype determined gender, see Sex column',
                   'SampleInfo',1,1,1,1,0,1,1,1, 'F=Female;M=Male;?=Unknown')
  
  # Add QC column descriptions
  for(c in colnames(combined_metrics)[which(grepl('_alignment$', colnames(combined_metrics)))]){
    info <- add_info(info, c, 'PicardMetrics after alignment, see https://broadinstitute.github.io/picard/picard-metric-definitions.html for details',
                           'QC', 1,1,1,1,1,1,1,1)
  }
  for(c in colnames(combined_metrics)[which(grepl('_genotyping$', colnames(combined_metrics)))]){
    info <- add_info(info, c, 'PicardMetrics after MarkDuplicates (only for genotyping from RNAseq), see https://broadinstitute.github.io/picard/picard-metric-definitions.html for details',
                           'QC', 0,0,0,0,1,0,0,0)
  }
  
  for(c in colnames(combined_metrics)[which(grepl('_R2$|_R1$', colnames(combined_metrics)))]){
    info <- add_info(info, c,'FastQC','QC', 0,1,1,1,1,1,1,1)
  }
  
  
  colnames(info) <- c('Column name', 'Description', 'Column type','encoding','Total_NotAvailable',
                      'AMP-AD', 'Braineac', 'Brainseq', 'CMC', 'ENA', 'GTEx', 'NABEC', 'TargetALS')
  
  columns_not_in <- colnames(combined_metrics)[!colnames(combined_metrics) %in% info$`Column name`]
  if(length(columns_not_in) > 0){
    print(columns_not_in)
    stop(paste0(length(columns_not_in), ' columns without description'))
  }
  # Make order in info same as in results
  info <- info[match(colnames(combined_metrics), info$`Column name`),]
  return(info)
}