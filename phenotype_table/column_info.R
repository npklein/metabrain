# add description of columns in dataframe
get_column_info <- function(phenotype_data){
  #number_NA <- na_count(phenotype_data)
  # Add info on what each of the columns mean, and check that there is info included for all columns in the output file
  add_info <- function(df, colname, description, coltype, encoding=''){
    df <- rbind(df, data.frame(colname,description, coltype, encoding))
      return(df)
  }
  info <- data.frame(col1=c(), col2=c())
  
  # Add sample info column descriptions
  info <- add_info(info, 'rnaseq_id','RNAseq identifier that matches FastQC ID','SampleInfo')
  info <- add_info(info, 'SampleFull',
                   'Sample name as used in the Compute5 pipeline we ran. Often in form of <individual ID>_<Sample/run_ID>, but not always',
                   'SampleInfo')
  info <- add_info(info, 'cohort','Name of the cohort. This is the larger cohort, e.g. cohort PsychEncode consists of smaller cohorts BipSeq, BrainVEX and others, and cohort ENA consists of many cohorts',
                   'SampleInfo')
  info <- add_info(info, 'MetaCohort','Name of the larger cohort, e.g. cohort PsychEncode, which consists out of multiple smaller cohorts. For smaller cohorts, see cohort column',
                   'SampleInfo')
  info <- add_info(info, 'Gender','Reported gender of the sample.',
                   'SampleInfo', 'F=Female;M=Male;')
  info <- add_info(info, 'educ','Education level. Highest grade or year of regular school as recorded during the baseline cognitive testing. Years of formal education was determined with the education question from the 1990 US Census',
                   'SampleInfo', 'Elementary 0 1 2 3 4 5 6 7 8; High School 9 10 11 12; College 13 14 15 16; Graduate\\Professional 17 18 19 20 21;')
  info <- add_info(info, 'AgeAtDiagnosis','Age of diagnosis','SampleInfo')
  info <- add_info(info, 'Age','Age at visit, possibly age at which sample is taken (not confirmed)','SampleInfo')
  info <- add_info(info, 'Weight_(pounds)','Weight in pound','SampleInfo')
  info <- add_info(info, 'Height_(Inches)','Height in inches','SampleInfo')
  info <- add_info(info, 'Brain_Weight_(in_grams)','Post mortem brain weight in grams','SampleInfo')
  info <- add_info(info, 'CauseDeath','Cause of death','SampleInfo')
  info <- add_info(info, 'AgeDeath','Age of death','SampleInfo')
  info <- add_info(info, 'library_layout','Library layout','SampleInfo')
  info <- add_info(info, 'libraryPrep','Library preparation method','SampleInfo')
  info <- add_info(info, 'RIN','RNA integrity number','SampleInfo')
  info <- add_info(info, 'PMI_(in_hours)','PMI in hours','SampleInfo')
  info <- add_info(info, 'pH','pH of postmortem brain','SampleInfo')
  info <- add_info(info, 'PMI_(in_hours)','Post Mortem Interval in hours','SampleInfo','SampleInfo')
  info <- add_info(info, 'cts_mmse30','The Mini Mental State Examination. A widely used, 30 item, standardized screening measure of dementia severity (see https://www.synapse.org/#!Synapse:syn3191090 for details)','SampleInfo')
  info <- add_info(info, 'braaksc','Braak Stage. This assessment is a semiquantitative measure of neurofibrillary tangles. Diagnosis includes algorithm and neuropathologists opinion.',
                   'SampleInfo',
                   '0 0;1 I;2 II;3 III;4 IV;5 V;6 VI;8 DK;9 Missing')
  info <- add_info(info, 'ceradsc','Assessment of neuritic plaques; Cross-sectional.  his assessment is a semiquantitative measure of neuritic plaques. Diagnosis includes algorithm and neuropathologists opinion',
                   'SampleInfo',
                   '1 Definite;2 Probable;3 Possible;4 No AD;9 Missing')
  info <- add_info(info, 'cogdx','Final Clinical Dx - Clinical Consensus Diagnosi.;Physicians overall cognitive diagnostic category','SampleInfo',
                   paste0('1 NCI, No cognitive impairment (No impaired domains);2 MCI, Mild cognitive impairment (One impaired domain) and NO other cause of CI;',
                          '3 MCI, Mild cognitive impairment (One impaired domain) AND another cause of CI','4 AD, Alzheimers disease and NO other cause of CI (NINCDS PROB AD)',
                          '5 AD, Alzheimers disease AND another cause of CI (NINCDS POSS AD)','6 Other dementia. Other primary cause of dementia'))
  info <- add_info(info, 'CDR','Clinical Dementia Rating','SampleInfo',
                   '0 no dementia; 0.5 Questionable dementia (MCI); 1 Mild dementia; 2 Moderate Dementia; 3 Sever dementia; 4 Profound dementia; 5 Terminal Dementia')
  info <- add_info(info, 'NP.1','Neuropathology Category as measured by CERAD','SampleInfo',
                   '(1 Normal; 2 Definite AD; 3 probable AD; 4 possible AD)')
  info <- add_info(info, 'PlaqueMean','Mean neocortical plaque density across 5 regions (# of plaques/mm2)','SampleInfo')
  info <- add_info(info, 'AgeOnset','Age of disease onset in years','SampleInfo')
  info <- add_info(info, 'Diagnosis','Disease diagnosis','SampleInfo')
  info <- add_info(info, 'OtherDiagnosis','Second, third etc diagnosis','SampleInfo')
  
  info <- add_info(info, 'FamilyHistory',' 	Family history of mental disorders','SampleInfo')
  info <- add_info(info, 'Lifetime Antipsychotics','Lifetime of antipsychotics in fluphenazine equivalents (mg)','SampleInfo')
  info <- add_info(info, 'History of neurological disorders','Stroke, epilepsy,neurodegenerative dz (Parkinsons dz, dementia, etc)','SampleInfo')
  for(f in c("Ethanol","Delta-9_THC_(Active)","Delta-9_Carboxy_THC_(Inactive)","11-Hydroxy_Delta-9_THC_(Active)","Cocaine","Benzoylecgonine","Cocaethylene","Nicotine","Cotinine","Codeine","Hydrocodone","Hydromorphone","Morphine","Fentanyl","Oxycodone","Oxymorphone","6-AM","Tramadol","Methadone","Benzodiazepines","PCP","Anticonvulsants","Antidepressants/_SSRIs","Antipsychotics","Antidepressants")){
    info <- add_info(info, f,'Toxicology test','SampleInfo')
  }
  
  info <- add_info(info, 'MedicalHistory (Neuro, Vascular)','Info on medical history','SampleInfo')
  info <- add_info(info, 'individualID','ID that corresponds to the individual.','SampleInfo')
  info <- add_info(info, 'genotype_id','ID that corresponds to the ID in the genotype file that we use.','SampleInfo')
  info <- add_info(info, 'SpecificBrainRegion','Specific brain region.','SampleInfo')
  info <- add_info(info, 'SpecificBrainRegion','Specific brain region.','SampleInfo')
  info <- add_info(info, 'apoe_genotype','apolipoprotein E (APOE)','SampleInfo','22 E2E2; 23 E2E3; 24; E2E4; 33; E3E3; 34 E3E4; 44 E4E4')
  info <- add_info(info, 'BrodmannArea','Brodmann area from which RNAseq sample is taken','SampleInfo')
  info <- add_info(info, 'Hemisphere','Hemisphere dissected','SampleInfo')
  info <- add_info(info, 'RNA_isolation_28S_18S_ratio','RNA isolation','SampleInfo')
  info <- add_info(info, 'RNA_isolation_TotalYield_ug','RNA isolation','SampleInfo')
  info <- add_info(info, 'RNA_isolation_260_280_ratio','RNA isolation','SampleInfo')
  info <- add_info(info, 'Total_DNA_ug','DNA isolation','SampleInfo')
  info <- add_info(info, 'DNA_isolation_Total_Yield_ug','DNA isolation','SampleInfo')
  info <- add_info(info, 'DNA_isolation_260_280_ratio','DNA isolation','SampleInfo')
  info <- add_info(info, 'DNA_isolation_260_230_ratio','DNA isolation','SampleInfo')
  info <- add_info(info, 'library_selection','RNA selection method','SampleInfo')
  info <- add_info(info, 'instrument_model','Sequencing model (should be RNA sequencing model, not always clear)','SampleInfo')
  info <- add_info(info, 'BroadBrainRegion','Larger brain region (e.g. Cortex, Cerebllum, Spincal cord)','SampleInfo')
  info <- add_info(info, 'Family_History_of_ALS/FTD?','Family history of ALS/FTD','SampleInfo')
  info <- add_info(info, 'Site_of_Motor_Onset','Site of motor onset','SampleInfo')
  info <- add_info(info, 'Site_of_Motor_Onset_Detail','Detailed site of motor onset','SampleInfo')
  info <- add_info(info, 'Disease_Duration_(Onset_to_Tracheostomy_or_Death)_in_Months','Detailed duration in months','SampleInfo')
  info <- add_info(info, 'Reported_Genomic_Mutations','Reported genomic mutation (such as c9orf repeat epxansion)','SampleInfo')
  info <- add_info(info, 'has_C9orf27_repeat_expansion','Yes/No has C9orf72 repeat expansion','SampleInfo')
  info <- add_info(info, 'has_ATXN2_repeat_expansion','Yes/No has ATXN2 repeat expansion','SampleInfo')
  info <- add_info(info, 'C9_repeat_size','Size of the C9orf72 repeat expansion','SampleInfo')
  info <- add_info(info, 'ATXN2_repeat_size','Size of the ATXN2 repeat expansion','SampleInfo')
  info <- add_info(info, 'Seizure_notes','Seizure notes','SampleInfo')
  info <- add_info(info, 'Seizures','Yes/No seizures','SampleInfo')
  info <- add_info(info, 'Comorbidities','Other comorbidities','SampleInfo')
  info <- add_info(info, 'RnaSequencingPlatform','The RNA sequencing platform. Beware: although I think they are all the rnaseq platform, some might be the dna sequencing plator','SampleInfo')
  info <- add_info(info, 'Modifier_of_ALS_Spectrum_MND_-_Family_History_of_ALS/FTD?','Modifier of ALS','SampleInfo')
  info <- add_info(info, 'Modifier_of_ALS_Spectrum_MND_-_FTD?','Modifier of ALS','SampleInfo')
  info <- add_info(info, 'Modifier_of_ALS_Spectrum_MND_-_Dementia?','Modifier of ALS','SampleInfo')
  info <- add_info(info, 'ERCC_Added','ERCC spike-in added','SampleInfo')
  info <- add_info(info, 'TissueState','State tissue was kept in','SampleInfo')
  info <- add_info(info, 'YearAutopsy','Year of the autopsy','SampleInfo')
  info <- add_info(info, 'MedicationsDeath','Medications taken at time of death','SampleInfo')
  info <- add_info(info, 'Smoker','Yes/No smoker','SampleInfo')
  info <- add_info(info, 'History_of_neurological_disorders','History of neurological disorders','SampleInfo')
  info <- add_info(info, 'EtohResults','Alcohol consumption','SampleInfo')
  info <- add_info(info, 'Detailed.Diagnosis','Extra information diagnosis','SampleInfo')
  info <- add_info(info, 'MedicalHistory_(Neuro,_Vascular)','Medical history','SampleInfo')
  info <- add_info(info, 'HealthConditions','Health condition','SampleInfo')
  info <- add_info(info, 'Medication_notes','Medication notes','SampleInfo')
  info <- add_info(info, 'Lifetime_Antipsychotics','Lifetime of antipsychotics in fluphenazine equivalents (mg)','SampleInfo')
  info <- add_info(info, 'IQ','IQ','SampleInfo')
  info <- add_info(info, 'IQ.notes','IQ notes','SampleInfo')
  info <- add_info(info, 'Pyschiatric.Medications','Yes/No History of psychriatic medications','SampleInfo')
  info <- add_info(info, 'FSIQ','Full Scale Intelligence Quotient','SampleInfo')
  
  
  info <- add_info(info, 'ADI.R.A..cut.off.10.','ADI-R A - score > 10 for autism','SampleInfo')
  info <- add_info(info, 'ADI.R.B..NV..cut.off.7.','ADI-R B - for nonverbal individuals - score > 7 for autism','SampleInfo')
  info <- add_info(info, 'ADI.R.D..cut.off.1.','ADI-R B - for verbal individuals - score > 1 for autism','SampleInfo')
  info <- add_info(info, 'ADI.R.B..V..cut.off.8.','ADI-R B - for verbal individuals - score > 8 for autism','SampleInfo')
  info <- add_info(info, 'ADI.R.C..cut.off.3.','ADI-R D- score > 3 for autism','SampleInfo')
  info <- add_info(info, 'Agonal.State','Cause of death categorized into agonal state as classified into','SampleInfo',
                   'S short; I intermediate; P prolonged; SZ related to seizures')
  
    # Add QC column descriptions
  #for(c in colnames(phenotype_data)[which(grepl('_alignment$', colnames(phenotype_data)))]){
  #  info <- add_info(info, c, 'PicardMetrics after alignment, see https://broadinstitute.github.io/picard/picard-metric-definitions.html for details',
  #                         'QC')
  #}
  #for(c in colnames(phenotype_data)[which(grepl('_genotyping$', colnames(phenotype_data)))]){
  #  info <- add_info(info, c, 'PicardMetrics after MarkDuplicates (only for genotyping from RNAseq), see https://broadinstitute.github.io/picard/picard-metric-definitions.html for details',
  #                         'QC')
  #}
  
  #for(c in colnames(phenotype_data)[which(grepl('_R2$|_R1$', colnames(phenotype_data)))]){
  #  info <- add_info(info, c,'FastQC','QC')
  #}
  
  
  #colnames(info) <- c('Column name', 'Description', 'Column type','encoding','Total_NotAvailable',
  #                    'AMP-AD', 'Braineac', 'Brainseq', 'CMC', 'ENA', 'GTEx', 'NABEC', 'TargetALS')
  phenotype_data <- phenotype_data[,colnames(phenotype_data) != 'column_to_keep']  
  columns_not_in <- colnames(phenotype_data)[!colnames(phenotype_data) %in% info$colname]
  if(length(columns_not_in) > 0){
    print(columns_not_in)
    stop(paste0(length(columns_not_in), ' columns without description'))
  }
  # Make order in info same as in results
  info <- info[match(colnames(phenotype_data), info$colname),]
  return(info)
}
