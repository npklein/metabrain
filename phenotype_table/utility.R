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
                metavar='character'),
    make_option(c('-g', '--genotype_IDs_dir'), type='character', default=getwd(), 
                help='Directory with files that contains genotype IDs', 
                metavar='character')
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
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
  # What somtimes happens with the merging of files is that there are multiple columns that should be the same
  # e.g. pH.x and pH.y (should be pH), or Daily Smoking and Smoking Per Day (should be Smoking).
  # This function takes 1+ columns in column_to_merge_into (e.g. pH.x and pH.y) and merges them into column_to_keep (e.g. pH)
  # Before it does this, it checks that the values in column_to_keep are NA. The reason for this is that if for example
  # a dataset has a column Smoking and a column Smoking Per Day, we don't want to overwrite Smoking with the Smoking Per Day value
  # This is a consequence of different cohorts having slightly different definitions which we would like to combine

  # Check if column_to_keep is NA, if so loop over columns in column_to_merge_into and try to itteratively fill it
  if(!is.data.table(dt)){
    stop('Error: input to change_rows_if_na has to be data.table')
  }
  # first check if not multiple columns are not NA 
  if(nrow(dt[length(column_to_merge_into)-rowSums(is.na(dt[,..column_to_merge_into]))>1,])){
    
    # multiple columns not NA, have to implement this first
    if(!overwrite){
      print(dt[length(column_to_merge_into)-rowSums(is.na(dt[,..column_to_merge_into]))>1,..column_to_merge_into])
      stop('multiple column not NA, implement this first')
    }
  }
  # pH column had trouble with converting implicitly when running on the cluster (works without this if statement on my laptop)
  # Because it will warn about type conversion. Therefore, set warning lower for this loop
  options(warn=0)
  for(column_name in column_to_merge_into){
    
    # Get the rows where column_to_keep == NA
    rows_na <- which(is.na(dt[,..column_name]) | dt[,..column_name] == '')
    # Then fill in with the value of the next column in the loop (might be NA, but doesn't matter). Do this for all columns to merge
    dt[rows_na, column_to_keep := dt[rows_na,][[column_name]]]
  }
  options(warn=2)
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
  dt[organism_part == 'Anterior cingulate cortex (BA24)',BrodmannArea := 'ba24']
  dt[body_site_annotations == 'Brain - Frontal Cortex (BA9)', BrodmannArea := 'ba9']
  dt[organism_part == 'Brodmann (1909) area 24', organism_part := NA ]
  dt[organism_part == 'Brodmann (1909) area 9', organism_part := NA ]
  dt[organism_part == 'Anterior cingulate cortex (BA24)',organism_part := NA]
  dt[body_site_annotations == 'Brain - Anterior cingulate cortex (BA24)',body_site_annotations := NA]
  dt[body_site_annotations == 'Brain - Frontal Cortex (BA9)', body_site_annotations := NA]
  
  
  t <- apply(dt, 1, function(r) any(r %in% c("Brain_-_Frontal_Cortex_(BA9)")))
  return(dt)
}
  
convert_to_broad_region <- function(dt){
  dt[,SpecificBrainRegion := gsub(' ','_',SpecificBrainRegion)]
  
  dt[SpecificBrainRegion %in% c('DLPFC','dorsolateral_prefrontal_cortex'), SpecificBrainRegion := 'Dorsolateral_Prefrontal_Cortex']
  dt[SpecificBrainRegion %in% c('temporal_cortex','TemporalCortex'), SpecificBrainRegion := 'Temporal_Cortex']
  dt[SpecificBrainRegion %in% c('amygdala'), SpecificBrainRegion := 'Amygdala']
  dt[SpecificBrainRegion %in% c('cerebellum'), SpecificBrainRegion := 'Cerebellum']
  dt[SpecificBrainRegion %in% c('Cortex_Frontal','frontal_cortex'), SpecificBrainRegion := 'Cerebral_Frontal_Cortex']
  dt[SpecificBrainRegion %in% c('HIPP'), SpecificBrainRegion := 'hippocampus']
  dt[SpecificBrainRegion %in% c('PUTM','Putamen_(basal_ganglia)'), SpecificBrainRegion := 'putamen']
  dt[SpecificBrainRegion %in% c('SNIG'), SpecificBrainRegion := 'substantia_nigra']
  dt[SpecificBrainRegion %in% c('sACC'), SpecificBrainRegion := 'subgenual_anterior_cingulate_cortex']
  dt[SpecificBrainRegion %in% c('Nucleus_accumbens_(basal_ganglia)','Nucleus_Accumbens'), SpecificBrainRegion := 'nucleus_accumbens']
  
  CX <- c('ba10','ba36','ba44','ba22','ba24', 'ba9',"cerebral_cortex","Cortex","Cortex_Occipital","Cortex_Sensory","Dorsolateral_Prefrontal_Cortex","Motor_Cortex",
          "Motor_Cortex_Lateral","Motor_Cortex_Medial","parietal_cortex","Prefrontal_Cortex","Temporal_Cortex","visual_cortex",'subgenual_anterior_cingulate_cortex',
          'Cerebral_Frontal_Cortex')
  
  CER <- c('Cerebellum','cerebellar_hemisphere')
  LIMBIC_SYSTEM <- c('hippocampus','hippocampus_proper','hypothalamus','Amygdala')
  BASAL_GANGLIA <- c('substantia_nigra','caudate_nucleus','putamen','PUTM_SNIG','nucleus_accumbens')
  SPINAL_CORD <- c("C1_segment_of_cervical_spinal_cord","Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Spinal_Cord_Thoracic" )
  dt$BroadBrainRegion <- ''
  dt[SpecificBrainRegion %in% CX, BroadBrainRegion := 'Cortex']
  dt[SpecificBrainRegion %in% CER, BroadBrainRegion := 'Cerebellum']
  dt[SpecificBrainRegion %in% LIMBIC_SYSTEM, BroadBrainRegion := 'Limbic system']
  dt[SpecificBrainRegion %in% BASAL_GANGLIA, BroadBrainRegion := 'Basal ganglia']
  dt[SpecificBrainRegion %in% SPINAL_CORD, BroadBrainRegion := 'Spinal cord']
  dt[SpecificBrainRegion=='Liver', BroadBrainRegion:='Liver']
  dt[grepl('Cell_Line', SpecificBrainRegion), BroadBrainRegion:='Cell line']
  dt[BroadBrainRegion=='',BroadBrainRegion:=NA]
  
  
  return(dt)
}
  
  
harmonize_and_clean_col_values <- function(merged_samplesheets){
    merged_samplesheets_fixed_columns <- copy(merged_samplesheets)
    ##### Remove some columns by name and all columns that have only NA values ####
    merged_samplesheets_fixed_columns[,c('Ethnicity','Race','race','RACE.inferred','Ancestry','Flowcell',
                                         'spanish','V47','V48','WGS_Source_Tissue_Type','organ','nucleicAcidSource',
                                         'BrainBank','BrainBank.x','BrainBank.y','FlowcellBatch',
                                         'pct_african','pct_south_asian','pct_east_asian','pct_european',
                                         'pct_americas','batch','V50','Source','Action','SnpArrayPlatform',
                                         'Genotyping_Chip_Type','Individual_Notes','Year_of_Autopsy','Source Name',
                                         'DNAmAge','Have you ever been a smoker?','AvgSpotLen','InsertSize',
                                         'is_tumor','submitted_subject_id','Site Sample Collected','individualIdSource',
                                         'isStranded','Material type','MND with FTD?','MND with Dementia?',
                                         'Expanded, >=30','cellType','Intermediate, 30-33','Revised El Escorial Criteria',
                                         'Cause of Death Related to ALS or MND?','PrimaryCauseOfDeath',
                                         'Other Family Notes'):=NULL]
    merged_samplesheets_fixed_columns <- merged_samplesheets_fixed_columns %>% replace_with_na(condition = ~.x == "")
    merged_samplesheets_fixed_columns <-   merged_samplesheets_fixed_columns[!merged_samplesheets_fixed_columns$rnaseq_id==T,]
    
    #### For some is reported Gender, others Sex, or other names. Merge into one column and then harmonize ####
    change_rows_if_na(merged_samplesheets_fixed_columns,'Gender',c("Sex","Sequenced Sex","Sex Call","Sex Genotype",
                                                                   "Genomic Sex","Sex (Sequenced)", "Sex.x","Sex.y",
                                                                   "SEX.inferred"),overwrite=T)
    merged_samplesheets_fixed_columns[Gender=='female' | Gender == 'Female', Gender := 'F']
    merged_samplesheets_fixed_columns[Gender=='male' | Gender == 'Male', Gender := 'M']
    merged_samplesheets_fixed_columns[Gender=='?', Gender := NA]
    ##### fix the tissue/brainregion ####
    change_rows_if_na(merged_samplesheets_fixed_columns,'BrodmannArea',c("Brodmann_Area","BroadmannArea","brodmannArea"))
    merged_samplesheets_fixed_columns <- fix_broadman_area(merged_samplesheets_fixed_columns)
    merged_samplesheets_fixed_columns$SpecificBrainRegion <- ''
    change_rows_if_na(merged_samplesheets_fixed_columns, 'SpecificBrainRegion',
                      c('BrainRegion','Brain_region','tissue','Tissue','organism_part','body_site','body_site_annotations','Sample Source','BrodmannArea'),
                      overwrite=T,
                      no_delete=c('BrodmannArea'))
    merged_samplesheets_fixed_columns <- convert_to_broad_region(merged_samplesheets_fixed_columns)
    
    ##### merge similar columns ####
    change_rows_if_na(merged_samplesheets_fixed_columns, 'SequencingPlatform',c('Instrument','platform','Platform','instrument_platform'),overwrite=T)
    setnames(merged_samplesheets_fixed_columns, old=c('SequencingPlatform'),new=c('RnaSequencingPlatform'))
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Diagnosis',c('Diagnosis.x','Diagnosis.y','Dx','Subject Group'))
    
    ##### Diagnosis ####
    merged_samplesheets_fixed_columns[Diagnosis=='Other Neurological Disorders' | Diagnosis == 'Other MND', Diagnosis := `Subject Group Subcategories`]
    merged_samplesheets_fixed_columns[Diagnosis=='SCZ' | Diagnosis == 'Schizo', Diagnosis := 'Schizophrenia']
    merged_samplesheets_fixed_columns[Diagnosis=='AD' | Diagnosis=='Alzheimer\'s Disease (AD)', Diagnosis := 'Alzheimer disease']
    merged_samplesheets_fixed_columns[Diagnosis=='AFF', Diagnosis := 'Affective Disorder']
    merged_samplesheets_fixed_columns[Diagnosis=='BP', Diagnosis := 'Bipolar']
    merged_samplesheets_fixed_columns[Diagnosis=='ASD', Diagnosis := 'Autism Spectrum Disorder']
    merged_samplesheets_fixed_columns[Diagnosis=='BP (not BP)' | Diagnosis==T | Diagnosis=='undetermined', Diagnosis := NA]
    merged_samplesheets_fixed_columns[Diagnosis=='PSP', Diagnosis := 'Progressive Supranuclear Palsy']
    merged_samplesheets_fixed_columns[Diagnosis=='ALS Spectrum MND, Other Neurological Disorders', Diagnosis := 'ALS Spectrum MND and Other Neurological Disorders']
    merged_samplesheets_fixed_columns[Diagnosis=='Control', Diagnosis := 'Non-Neurological Control']
    merged_samplesheets_fixed_columns[Diagnosis=='ALS Spectrum MND and Other Neurological Disorders', Diagnosis := 'ALS Spectrum MND']
    merged_samplesheets_fixed_columns[Diagnosis=='Pre-fALS, Other Neurological Disorders', Diagnosis := 'Pre-fALS']
    
    # There is sometimes overlap between subcategories and subcategory, get those out
    setnames(merged_samplesheets_fixed_columns, old = c('Subject Group Subcategories'), new = c('OtherDiagnosis'))
    merged_samplesheets_fixed_columns[!is.na(OtherDiagnosis) & !(grepl('Classical/Typical ALS',OtherDiagnosis) & grepl('Classical/Typical ALS',`Subject Group Subcategory`)), 
                                      OtherDiagnosis := paste0(`Subject Group Subcategory`,'; ',OtherDiagnosis)]
    merged_samplesheets_fixed_columns[is.na(OtherDiagnosis), OtherDiagnosis := `Subject Group Subcategory` ]
    
    merged_samplesheets_fixed_columns[OtherDiagnosis=='Alzheimer\'s Disease (AD); Alzheimer\'s Disease (AD)', OtherDiagnosis := 'Alzheimer disease']
    merged_samplesheets_fixed_columns[OtherDiagnosis=='; Non-Neurological Control', OtherDiagnosis := NA]
    
    merged_samplesheets_fixed_columns[OtherDiagnosis=='Classical/Typical ALS, Alzheimer\'s Disease (AD), Chronic Traumatic Encephalopathy (CTE); Classical/Typical ALS ', 
                                      OtherDiagnosis := 'Classical/Typical ALS, Alzheimer\'s Disease (AD), Chronic Traumatic Encephalopathy (CTE); Classical/Typical ALS ']
    merged_samplesheets_fixed_columns[,`Subject Group Subcategory` := NULL]
    ######
    
    change_rows_if_na(merged_samplesheets_fixed_columns, 'libraryPrep',c('LibraryPrep','LibraryKit','Prep'),overwrite=T)
    merged_samplesheets_fixed_columns[grepl('Not Available|Unknown|?', pH), pH := NA]
    merged_samplesheets_fixed_columns[,pH := as.numeric(pH)]

    change_rows_if_na(merged_samplesheets_fixed_columns, 'pH',c('pH.x','pH.y'),overwrite=T)
    merged_samplesheets_fixed_columns[grepl('?', `PMI_(in_hours)`), `PMI_(in_hours)` := NA]
    merged_samplesheets_fixed_columns[,`PMI_(in_hours)` := as.numeric(`PMI_(in_hours)`)]
    merged_samplesheets_fixed_columns[grepl('Not Applicable|Unknown|?', `Post Mortem Interval in Hours`), `Post Mortem Interval in Hours` := NA]
    merged_samplesheets_fixed_columns[,`Post Mortem Interval in Hours` := as.numeric(`Post Mortem Interval in Hours`)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'PMI_(in_hours)',c('PMI.x','PMI.y','PMI','Post Mortem Interval in Hours'),overwrite=T)
    
    merged_samplesheets_fixed_columns[grepl('Not Available|Unknown|?', RIN), RIN := NA]
    merged_samplesheets_fixed_columns[, RIN := as.numeric(RIN)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'RIN',c('RIN.x','RIN.y'), overwrite=T)
    
    #### library layout ####
    change_rows_if_na(merged_samplesheets_fixed_columns, 'library_layout',c('runType','RunType'))
    merged_samplesheets_fixed_columns[library_layout=='PAIRED' | library_layout == 'pairedEnd', library_layout := 'paired-end']
    merged_samplesheets_fixed_columns[library_layout=='SINGLE',library_layout := 'single-end']
    
    
    
    ##### death #####
    change_rows_if_na(merged_samplesheets_fixed_columns, 'AgeDeath',c('age_death','AOD'))
    change_rows_if_na(merged_samplesheets_fixed_columns, 'CauseDeath',c('MannerOfDeath','Cause of Death','Manner_Of_Death','Cause_of_Death','CODE'))
    merged_samplesheets_fixed_columns[grepl('UNDETERMINED|Unknown|UNKOWN|Yes', CauseDeath), CauseDeath := NA ]
    merged_samplesheets_fixed_columns[CauseDeath == 'ASPHYXIA', CauseDeath := 'Asphyxia' ]
    merged_samplesheets_fixed_columns[CauseDeath == 'ASTHMA', CauseDeath := 'Asthma' ]
    merged_samplesheets_fixed_columns[CauseDeath == 'CARDIAC', CauseDeath := 'Cardiac' ]
    merged_samplesheets_fixed_columns[CauseDeath == 'Complications fro mALS' | CauseDeath == 'Complications from ALS', CauseDeath := 'Complications of ALS' ]
    
    ##### brain weight and other physical phenotypes ####
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Brain_Weight_(in_grams)',c('Brain Weight (gram)','BrainWeight'))
    merged_samplesheets_fixed_columns[grepl('N/A', Height), Height := NA]
    merged_samplesheets_fixed_columns[,Height:=as.numeric(Height)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Height_(Inches)',c('Height (inches)','Height'))
    merged_samplesheets_fixed_columns[grepl('N/A', Weight), Weight := NA]
    merged_samplesheets_fixed_columns[,Weight:=as.numeric(Weight)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Weight_(pounds)',c('Weight','Weight (lbs)'))    
  
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Age',c('age_at_visit_max'))
    merged_samplesheets_fixed_columns$AgeAtDiagnosis <- ''
    merged_samplesheets_fixed_columns[,`Age of Diagnosis`:=as.character(`Age of Diagnosis`)]
    merged_samplesheets_fixed_columns[,`Age at Diagnosis`:=as.character(`Age at Diagnosis`)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'AgeAtDiagnosis',c('Age at Diagnosis','Age of Diagnosis',
                                                                            'age_first_ad_dx'), overwrite=T)
    
    change_rows_if_na(merged_samplesheets_fixed_columns, 'cts_mmse30_first_ad_dx',c('cts_mmse30_lv'), overwrite=T)
    setnames(merged_samplesheets_fixed_columns, old=c('cts_mmse30_first_ad_dx'),new=c('cts_mmse30'))
    
    merged_samplesheets_fixed_columns[,braaksc := as.integer(braaksc)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'braaksc',c('bbscore'))
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Smoker', c('SmokingEither','Are you an active smoker?'))
    merged_samplesheets_fixed_columns[grepl('nknown',merged_samplesheets_fixed_columns$Smoker), Smoker := NA]
    merged_samplesheets_fixed_columns[Smoker=='no', Smoker := 'No']
    merged_samplesheets_fixed_columns[Smoker=='yes', Smoker := 'Yes']

    #### age onset ####
    merged_samplesheets_fixed_columns[grepl('Not Applicable|Unknown', `Age at Symptom Onset`), `Age at Symptom Onset` := NA]
    merged_samplesheets_fixed_columns[grepl('Not Applicable|Unknown', `Age of Symptom Onset`), `Age of Symptom Onset` := NA]
    merged_samplesheets_fixed_columns[, `Age of Symptom Onset` := as.integer(`Age of Symptom Onset`)]
    merged_samplesheets_fixed_columns[, `Age at Symptom Onset` := as.integer(`Age at Symptom Onset`)]

    change_rows_if_na(merged_samplesheets_fixed_columns, 'AgeOnset', c('Age of Symptom Onset','AgeOnsetMania','Age at Symptom Onset','AgeOnsetSchizo'),overwrite=T)
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Hemisphere', c('hemisphere'))
    
    merged_samplesheets_fixed_columns[`Family History of ALS/FTD?`=='Not Applicable' | `Family History of ALS/FTD?` == 'Unknown', `Family History of ALS/FTD?` := NA ]

    change_rows_if_na(merged_samplesheets_fixed_columns, 'Site of Motor Onset', c('Site of Motor'))
    merged_samplesheets_fixed_columns[, DurationIllness := DurationIllness * 12]
    merged_samplesheets_fixed_columns[`Disease Duration (Onset to Tracheostomy or Death) in Months` == 'Unknown', `Disease Duration (Onset to Tracheostomy or Death) in Months` := NA]
    merged_samplesheets_fixed_columns[`Disease Duration (Onset to Tracheostomy or Death) in Months` == 'Not Applicable', `Disease Duration (Onset to Tracheostomy or Death) in Months` := NA]
    merged_samplesheets_fixed_columns[,`Disease Duration (Onset to Tracheostomy or Death) in Months` := gsub('~','', `Disease Duration (Onset to Tracheostomy or Death) in Months`)]
    merged_samplesheets_fixed_columns[,`Disease Duration (Onset to Tracheostomy or Death) in Months` := as.numeric(`Disease Duration (Onset to Tracheostomy or Death) in Months`)]
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Disease Duration (Onset to Tracheostomy or Death) in Months', c('DurationIllness'))
    
    change_rows_if_na(merged_samplesheets_fixed_columns, 'Reported Genomic Mutations', c('Reported Genomic Mutations (from Site, not associated with sequencing data from NYGC) '))
    
    merged_samplesheets_fixed_columns[`C9 repeat size` == 'N/A' |`C9 repeat size` == '', `C9 repeat size` := NA ]
    merged_samplesheets_fixed_columns[`C9 Expansion per IGM (TALS only)` == 'Unknown', `C9 Expansion per IGM (TALS only)` := NA ]
    merged_samplesheets_fixed_columns[`C9orf72 Repeat Expansion (CUMC, Sep 2018)` == 'Unknown', `C9orf72 Repeat Expansion (CUMC, Sep 2018)` := NA ]
    
    merged_samplesheets_fixed_columns[`ATXN2 repeat size` == 'N/A' |`ATXN2 repeat size` == '', `ATXN2 repeat size` := NA ]
    merged_samplesheets_fixed_columns[`ATXN2 Expansion per IGM (TALS only)` == 'Unknown', `ATXN2 Expansion per IGM (TALS only)` := NA ]
    merged_samplesheets_fixed_columns[`ATXN2 Repeat Expansion? (CUMC, Sep2018)` == 'Unknown', `ATXN2 Repeat Expansion? (CUMC, Sep2018)` := NA ]
    
    change_rows_if_na(merged_samplesheets_fixed_columns, 'C9 Expansion per IGM (TALS only)', c('C9orf72 Repeat Expansion (CUMC, Sep 2018)'))
    change_rows_if_na(merged_samplesheets_fixed_columns, 'ATXN2 Expansion per IGM (TALS only)', c('ATXN2 Repeat Expansion? (CUMC, Sep2018)'))
    
    setnames(merged_samplesheets_fixed_columns, old = c('C9 Expansion per IGM (TALS only)','ATXN2 Expansion per IGM (TALS only)'), 
             new = c('has_C9orf27_repeat_expansion','has_ATXN2_repeat_expansion'))
    
    merged_samplesheets_fixed_columns[Seizure_notes == 'Epiepsy' | Seizure_notes == 'Epliepsy' | Seizure_notes=='Epilepsy - Seizure disorder, temporal lobe', Seizure_notes := 'Epilepsy' ]
    merged_samplesheets_fixed_columns[Seizure_notes=='no information' | Seizure_notes == 'no records' | Seizure_notes=='',Seizure_notes:=NA]
    merged_samplesheets_fixed_columns[Seizure_notes=='noted no significant history' | Seizure_notes =='noted to have no significant medical history', Seizure_notes:=  'No notable medical history']
    
    #### comorbidities ####
    merged_samplesheets_fixed_columns[Comorbidities=='Unknown', Comorbidities:=NA]
    merged_samplesheets_fixed_columns[Comorbidity.notes..other.than.seizures.=='',Comorbidity.notes..other.than.seizures.:=NA]
    merged_samplesheets_fixed_columns[is.na(Comorbidities), Comorbidities := Comorbidity.notes..other.than.seizures.]
    merged_samplesheets_fixed_columns[!is.na(Comorbidities) & !is.na(Comorbidity.notes..other.than.seizures.), Comorbidities := paste0(Comorbidities,'; ',Comorbidity.notes..other.than.seizures.)]
    merged_samplesheets_fixed_columns[,Comorbidity.notes..other.than.seizures.:=NULL]
    
    
    #### modifiers #####
    merged_samplesheets_fixed_columns[`Modifier of ALS Spectrum MND - Family History of ALS/FTD?`=='Unknown' | `Modifier of ALS Spectrum MND - Family History of ALS/FTD?` == 'Not Applicable', `Modifier of ALS Spectrum MND - Family History of ALS/FTD?` := NA]
    merged_samplesheets_fixed_columns[`Modifier of ALS Spectrum MND - FTD?`=='Unknown' | `Modifier of ALS Spectrum MND - FTD?` == 'Not Applicable', `Modifier of ALS Spectrum MND - FTD?` := NA]
    merged_samplesheets_fixed_columns[`Modifier of ALS Spectrum MND - Dementia?`=='Unknown' | `Modifier of ALS Spectrum MND - Dementia?` == 'Not Applicable', `Modifier of ALS Spectrum MND - Dementia?` := NA]

    merged_samplesheets_fixed_columns[IQ.notes=='', IQ.notes:=NA]
    
    table(merged_samplesheets_fixed_columns$ADI.R.A..cut.off.10.)
                                                      
    return(merged_samplesheets_fixed_columns)
  }
