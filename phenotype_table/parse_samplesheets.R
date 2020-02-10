# Parse samplesheets of the different Brain datasets
library(data.table)
source('utility.R')
library(naniar)
library(stringr)

parse_samplesheet_files <- function(opt){
  AMP_AD <- parse_amp_ad(opt)
  Braineac <- parse_Braineac(opt)
  Brainseq <- parse_Brainseq(opt)
  CMC <- parse_CMC(opt)
  GTEx <- parse_GTEx(opt)
  NABEC <- parse_NABEC(opt)
  TargetALS <- parse_TargetALS(opt)
  ENA <- parse_ENA(opt)
  BrainGVEx <- parse_BrainGVEx(opt)
  UCLA_ASD <- parse_UCLA_ASD(opt)
  BipSeq <- parse_BipSeq(opt)
  CMC_HBCC <- parse_CMC_HBCC(opt)
  merged_samplesheets <- data.frame()
  for(df in list(AMP_AD, Braineac, Brainseq, CMC, GTEx,
                 NABEC, TargetALS, ENA, BrainGVEx,UCLA_ASD,
                 BipSeq,CMC_HBCC)){
      merged_samplesheets <- rbind(merged_samplesheets, df, fill=T)
  }
  # remove some columns that are easier to remove now than for specific cohorts
  merged_samplesheets <- harmonize_and_clean_col_values(merged_samplesheets)
  return(merged_samplesheets)
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
  setnames(ROSMAP, old=c('wgs_id','study', 'msex','pmi'), new=c('genotype_id','cohort','Gender','PMI'))
  
  # age_at_visit is empty, remove
  ROSMAP <- ROSMAP[,c('projid','age_at_visit'):=NULL]
  ROSMAP <- ROSMAP[,which(unlist(lapply(ROSMAP, function(x)!all(is.na(x))))),with=F]

  # Change sex code to M/F to match other cohorts
  ROSMAP[, Gender := as.character(Gender)][Gender == "1", Gender := "M"]
  ROSMAP[Gender == "0", Gender := "F"]
  
  # Change empty values to NA for easier processing later
  ROSMAP[which(ROSMAP$age_first_ad_dx==""),]$age_first_ad_dx <- NA
  ROSMAP <- ROSMAP[!duplicated(ROSMAP$rnaseq_id)]
  ROSMAP$Tissue <- 'DLPFC'
  return(ROSMAP)
}

parse_Mayo <- function(opt){
  Mayo_wgs <- fread(paste0(opt$metadataDir,'/AMP_AD/AMP-AD_Mayo_WGS_covariates.csv'))
  Mayo_wgs <- Mayo_wgs[,c('WGS_Participant_ID','Diagnosis','Sex','AgeAtDeath','ApoE','PMI','WGS_Source_Tissue_Type',
                          'RNAseq_CER_SampleID','RNAseq_TCX_SampleID')]
  
  Mayo_RNAseq_CBE <- fread(paste0(opt$metadataDir,'/AMP_AD/MayoRNAseq_RNAseq_CBE_covariates.csv'))
  Mayo_RNAseq_TCX <- fread(paste0(opt$metadataDir,'/AMP_AD/MayoRNAseq_RNAseq_TCX_covariates.csv'))
  
  setnames(Mayo_RNAseq_TCX, old = c('ID','FLOWCELL'), new = c('rnaseq_id', 'Flowcell'))
  setnames(Mayo_RNAseq_CBE, old = c('Sex'), new = c('Gender'))
  setnames(Mayo_wgs, old = c('Sex'), new = c('Gender'))

  colnames(Mayo_RNAseq_CBE)[colnames(Mayo_RNAseq_CBE)=='SampleID'] <- 'rnaseq_id'
  if(sum(Mayo_RNAseq_CBE$SampleID %in% Mayo_RNAseq_TCX$SampleID) > 0){ stop ('Should only be different SampleID in CBE and TCX files')}

  # The WGS data has 1 column for CBE samples, and one column for TCX smaple IDs, first merge one, then the other
  colnames(Mayo_wgs)[colnames(Mayo_wgs)=='RNAseq_CER_SampleID'] <- 'rnaseq_id'
  Mayo_CER <- merge(Mayo_RNAseq_CBE, Mayo_wgs, by=c('rnaseq_id', 'Diagnosis','AgeAtDeath','ApoE','PMI','Gender'),
                    all.x=T)
  same_rows(Mayo_CER, list(Mayo_RNAseq_CBE))
  Mayo_CER$RNAseq_TCX_SampleID <- NULL
  # for some reason for some samples it is NA, fix this
  Mayo_wgs$rnaseq_id <- NULL
  colnames(Mayo_wgs)[colnames(Mayo_wgs)=='RNAseq_TCX_SampleID'] <- 'rnaseq_id'
  Mayo_TCX <- merge(Mayo_RNAseq_TCX, Mayo_wgs, by=c('rnaseq_id', 'Diagnosis','AgeAtDeath','ApoE','PMI','Gender'),
                   all.x=T)
  same_rows(Mayo_TCX, list(Mayo_RNAseq_TCX))
  Mayo_CER$Tissue <- 'Cerebellum'
  Mayo_TCX$Tissue <- 'Cortex'
  Mayo <- rbind(Mayo_CER, Mayo_TCX)
  if(nrow(Mayo) != nrow(Mayo_RNAseq_CBE)+nrow(Mayo_RNAseq_TCX)){stop("Mayo should have equal number of rows as Mayo_RNAseq_CBE+Mayo_RNAseq_TCX")}
  
  setnames(Mayo, old = c('WGS_Participant_ID','ApoE','AgeAtDeath'), new = c('genotype_id', 'apoe_genotype','age_death'))
  Mayo[Mayo$age_death=='90_or_above',]$age_death <- '90+'
  Mayo$cohort <- 'Mayo'
  Mayo[, genotype_id := as.character(genotype_id)]
  Mayo[is.na(genotype_id),genotype_id := gsub('_CER|_TCX','',rnaseq_id)]
  
  return(Mayo)
}

parse_MSBB <- function(opt){
  MSBB_genotype_id_map <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_final.txt'),header=F)
  colnames(MSBB_genotype_id_map) <- c('rnaseq_id','genotype_id')

  MSBB_covariates <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_WES_covariatestsv'))
  colnames(MSBB_covariates)[colnames(MSBB_covariates)=='Region'] <- 'BrodmannArea'
  MSBB_covariates$barcode <- NULL
  MSBB_update <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_RNAseq_covariates_November2018Update.csv'))
  MSBB_update <- MSBB_update[MSBB_update$fileType=='fastq',]
  MSBB_update <- MSBB_update[,c('synapseId','fileName','barcode','fileType'):= NULL]
  
  if(nrow(MSBB_update) != length(unique(MSBB_update$sampleIdentifier))) { stop("should only be unique samples")}
  
  MSBB_covariates <- rbind.fill(MSBB_update, MSBB_covariates)
  MSBB_covariates$BrodmannArea <- gsub('-','',MSBB_covariates$BrodmannArea)
  
  MSBB_covariates$genotype_id <- MSBB_covariates$individualIdentifier.inferred
  MSBB_covariates[is.na(MSBB_covariates$individualIdentifier.inferred),]$genotype_id <- MSBB_covariates[is.na(MSBB_covariates$individualIdentifier.inferred),]$individualIdentifier
  MSBB_covariates[MSBB_covariates$genotype_id=='.',]$genotype_id <- MSBB_covariates[MSBB_covariates$genotype_id=='.',]$individualIdentifier
  
  MSBB_clinical <- fread(paste0(opt$metadataDir,'/AMP_AD/MSBB_clinical.csv'))
    
  MSBB <- merge(MSBB_covariates, MSBB_clinical, by.x='genotype_id', by.y='individualIdentifier',all.x=T)
  
  
  same_rows(MSBB, list(MSBB_covariates))
  colnames(MSBB)[colnames(MSBB)=='sampleIdentifier'] <- 'rnaseq_id'
  colnames(MSBB)[colnames(MSBB)=='RACE'] <- 'race'
  colnames(MSBB)[colnames(MSBB)=='SEX'] <- 'Gender'
  MSBB$cohort <- 'MSBB'
  
  # match race to ROSMAP race code (below from their column description)
  # MSBB:   A=Asian, B=Black, H=Hispanic, W=White.
  # ROSMAP: 1= White, 2=Black, Negro, African-American,3 Native American, Indian,4 Eskimo,5 Aleut,6 Asian or Pacific Island,8 REFUSAL,9 DON'T KNOW
  # And add 10 for Hispanic
  MSBB[MSBB$race=='A',]$race <- 6
  MSBB[MSBB$race=='B',]$race <- 2
  MSBB[MSBB$race=='H',]$race <- 10
  MSBB[MSBB$race=='W',]$race <- 1

  MSBB[MSBB$RACE.inferred=='A',]$RACE.inferred <- 6
  MSBB[MSBB$RACE.inferred=='B',]$RACE.inferred <- 2
  MSBB[MSBB$RACE.inferred=='H',]$RACE.inferred <- 10
  MSBB[MSBB$RACE.inferred=='W',]$RACE.inferred <- 1
  MSBB[MSBB$RACE.inferred=='.',]$RACE.inferred <- 9
  
  # make the apoe genotype the same as in the other cohorts
  MSBB$apoe_genotype <- NA
  MSBB[!is.na(MSBB$Apo1),]$apoe_genotype <- paste0(MSBB[!is.na(MSBB$Apo1),]$Apo1,MSBB[!is.na(MSBB$Apo1),]$Apo2)
  MSBB$Apo1 <- NULL
  MSBB$Apo2 <- NULL
  
  # Remove Mapped, TotalReads and rRNA rate because that should be in the Picard QC files that we also have
  MSBB$rRNA.rate <- NULL
  MSBB$TotalReads <- NULL
  MSBB$Mapped <- NULL
  
  MSBB$genotype_id <- NULL
  MSBB <- merge(MSBB, MSBB_genotype_id_map, by='rnaseq_id',all.x=T)
  MSBB$individualIdentifier <- MSBB$genotype_id
  MSBB$individualIdentifier.inferred <- NULL
  return(MSBB)
}

parse_amp_ad <- function(opt){
  Mayo <- parse_Mayo(opt)
  ROSMAP <- parse_ROSMAP(opt)
  MSBB <- parse_MSBB(opt)

  AMP_AD <- rbind.fill(Mayo, ROSMAP)
  AMP_AD <- rbind.fill(AMP_AD, MSBB)

  AMP_AD$SampleFull <- AMP_AD$rnaseq_id
  AMP_AD$MetaCohort <- 'AMP-AD'
  colnames(AMP_AD)[colnames(AMP_AD)=='PMI'] <- 'PMI_(in_hours)'
  AMP_AD$individualIdentifier <- NULL
  AMP_AD$individualIdentifier.inferred <- NULL
  
  return(AMP_AD)
}

parse_Braineac <- function(opt){
  Braineac <- fread(paste0(opt$metadataDir,'/Braineac/SampleInfoBiogen.csv'))
  setnames(Braineac, old=c('Sample_ID','SD No','Brain Bank','RIN-est','Region_simplified','PMI (h)'),
           new=c('rnaseq_id','genotype_id','cohort','RIN','Brain_region','PMI_(in_hours)'))
  
  Braineac$genotype_id <- gsub('/','_', Braineac$genotype_id)
  Braineac$MetaCohort <- 'Braineac'
  Braineac$SampleFull <- paste0(Braineac$genotype_id,'_',Braineac$rnaseq_id)
  Braineac$Individual_ID <- Braineac$genotype_id
  Braineac[rnaseq_id == 'A653-1341', rnaseq_id:='A653_1341']
  return(Braineac)
}

parse_Brainseq <- function(opt){
  Brainseq_phenotypes <- fread(paste0(opt$metadataDir,'/Brainseq/phenotypeFile_LIBD_szControl.csv'))
  Brainseq_phenotypes <- Brainseq_phenotypes[, .SD, .SDcols = ! names(Brainseq_phenotypes) %like% "snpPC"]
  
  # WGBS is methylation data, not needed (yet?)
  #Brainseq_WGBS_metadata <- fread(paste0(opt$metadataDir,'/Brainseq/LIBD_R21MH102791_LIBD-WGBS_Metadata.csv'))
  #Brainseq_WGBS_metadata <-   Brainseq_WGBS_metadata[,c('individualID', 'Cell.Type', 'Proportion.of.Neurons','yield.nuclei/mg.tissue',
  #                                                      'pg.DNA/nuclei.input','260/280.DNA','Data.Yield')]
  
  Brainseq_snpGenotype_metadata <- fread(paste0(opt$metadataDir,'/Brainseq/LIBD_R21MH102791_LIBD-WGBS_Metadata_SNPgenotypes_August2016Release.csv'))
  Brainseq_snpGenotype_metadata <- Brainseq_snpGenotype_metadata[,c('individualID','platform')]
  Brainseq_clinical_metadata <- fread(paste0(opt$metadataDir,'/Brainseq/LIBD_R21MH102791_LIBD-WGBS_Metadata_Clinical_August2016Release.csv'))
  Brainseq_clinical_metadata <-   Brainseq_clinical_metadata[,which(unlist(lapply(Brainseq_clinical_metadata, function(x)!all(is.na(x))))),with=F]
  
  Brainseq_clinical_metadata <-   Brainseq_clinical_metadata[,c('grant', 'study', 'Source','PMI'):=NULL]

  Brainseq <- merge(Brainseq_phenotypes, Brainseq_snpGenotype_metadata, by.x='RNum',by.y='individualID', all.x=T)
  
  Brainseq$cohort <- 'Brainseq'
  Brainseq$MetaCohort <- 'Brainseq'
  Brainseq$SampleFull <- paste0(Brainseq$individual_ID,'_',Brainseq$rnaseq_id)
  same_rows(Brainseq, list(Brainseq_phenotypes))
  
  Brainseq$rnaseq_id <- paste0(Brainseq$BrNum,'_',Brainseq$RNum)
  Brainseq$RNum <- NULL
  setnames(Brainseq, old=c('BrNum','Region','Sex','PMI'),
           new=c('genotype_id','Brain_region','Gender','PMI_(in_hours)'))
  return(Brainseq)
}

parse_CMC <- function(opt){
  CMC_clinical <- fread(paste0(opt$metadataDir,'/CMC/CMC-CMC_HBCC_clinical_.csv'))
  colnames(CMC_clinical) <- gsub(' ','_',colnames(CMC_clinical))
  
  # Remove CMC_HBCC samples, they are in separate file
  CMC_clinical <- CMC_clinical[!grepl('CMC_HBCC',CMC_clinical$Individual_ID),]
  
  CMC_rnaseq <- fread(paste0(opt$metadataDir,'/CMC/CMC_Human_DLPFC_rnaSeq.csv'))
  
  # This table does not seem to add much
  #CMC_rnaseq_synapse_table <- fread(paste0(opt$metadataDir,'/CMC/CMC_HBCC_RNA_table.txt'))
  #CMC_rnaseq_synapse_table <-   CMC_rnaseq_synapse_table[,which(unlist(lapply(CMC_rnaseq_synapse_table, function(x)!all(is.na(x))))),with=F]
  #CMC_rnaseq_synapse_table[,c('id','name','consortium','species','study','parentId','Contributor',
  #                                  'dataFileHandleId','currentVersion','createdBy','Capstone_1','assay','organ',
  #                                  'Capstone_4','nucleicAcidSource','fileFormat','isStranded','readLength','runType',
  #                            'RIN','notes','hemisphere','specimenID','tissue','BrodmannArea','platform') := NULL]
  #setnames(CMC_rnaseq_synapse_table, old=('individualID'), new=c('Individual_ID'))
  CMC_snp <- fread(paste0(opt$metadataDir,'/CMC/CMC_Human_SNP.csv'))
  CMC_dna_metadata <- fread(paste0(opt$metadataDir,'/CMC/CMC_MSSM-Penn-Pitt_DLPFC_DNA-metaData.csv'))
  CMC_rna_metadata <- fread(paste0(opt$metadataDir,'/CMC/CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv'))
  
  colnames(CMC_rna_metadata) <- gsub('DLPFC_','', colnames(CMC_rna_metadata))
  
  setnames(CMC_rna_metadata, old=c('RNA_Sequencing_Sample_ID','RNA_isolation_RIN','RNA_Sequencing_Ribozero_Batch',
                                   'RNA_Sequencing_Library_Batch','RNA_Sequencing_Flowcell_Batch'), 
                            new=c('rnaseq_id','RIN','Ribozero_Batch','Library_Batch','Flowcell_Batch'))
  
  # the following columsns are all based on mapping, which we redo. Therefore, these columns are removed. Also remove somef lowell columns
  CMC_rna_metadata <- CMC_rna_metadata[,c('RNA_Sequencing_Percent_Aligned','RNA_Sequencing_Total_Reads','RNA_Sequencing_rRNA_Rate',
                                          'RNA_Sequencing_Expression_Profiling_Efficiency','RNA_Sequencing_Transcripts_Detected',
                                          'RNA_Sequencing_Genes_Detected','RNA_Sequencing_Intergenic_Rate','RNA_Sequencing_Intronic_Rate',
                                          'RNA_Sequencing_Exonic_Rate','RNA_Sequencing_Mapped_Reads','RNA_Sequencing_Flowcell_Lane_B',
                                          'RNA_Sequencing_Flowcell_Lane_A','RNA_Sequencing_Flowcell_Name') := NULL]

  CMC_rna_metadata[, RNA_isolation_PrepOperator := as.character(RNA_isolation_PrepOperator)][RNA_isolation_PrepOperator == 1, RNA_isolation_PrepOperator := 'Jessica Johnson']
  CMC_rna_metadata[RNA_isolation_PrepOperator == 2, RNA_isolation_PrepOperator := 'Silvia De rubeis']
  setnames(CMC_rnaseq, old=c('Sample_RNA_ID', 'Flowcell_ID','Brain_Region', 'Operator',
                             'RNA_Prep_Date', '28S/18S','RNA_Prep_Operator','Total_RNA_Yield','260/280',
                             'RNA_Isolation_Batch'), 
                       new=c('rnaseq_id','Flowcell', 'Brain_region','Dissection_Operator',
                             'RNA_isolation_PrepDate','RNA_isolation_28S_18S_ratio','RNA_isolation_PrepOperator',
                             'RNA_isolation_TotalYield_ug','RNA_isolation_260_280_ratio','RNA_isolation_Batch'))
  
  # Remove a bunch of flowcell/lane information as it clutters the covariate/phenotype output. Should be able to see flowcell batch effects from the Flowcell column, if necesarry
  CMC_rnaseq <- CMC_rnaseq[,c('Expression_Profiling_Efficiency','Intergenic_Rate','Intronic_Rate','Intragenic_Rate','Percent_Aligned','Transcripts_Detected',
                              'Genes_Detected','Mapped_Reads','Total_Reads','Library_Barcode','Flowcell_(Given_to_Core)','Flowcell_Name','Flowcell_Name_2',
                              'Flowcell_Lane_A','Flowcell_Lane_B','Flowcell_2') := NULL]
  
  # only a few that are in metadata but not in rnaseq, rbind those to rnaseq
  CMC_rna_metadata_not_in_rnaseq <- CMC_rna_metadata[which(!CMC_rna_metadata$rnaseq_id %in% CMC_rnaseq$rnaseq_id),]
  CMC_rna_metadata_not_in_rnaseq$Flowcell <- NA
  CMC_rna_metadata_not_in_rnaseq$Study <- 'CMC'
  CMC_rna_metadata_not_in_rnaseq$Brodmann_Area <- NA
  CMC_rna_metadata_not_in_rnaseq$Brain_region <- NA
  CMC_rnaseq <- rbind(CMC_rnaseq , CMC_rna_metadata_not_in_rnaseq)
  
  # CMC_rnaseq_synapse_table doesn't add anything new from the other metadata files as far as I can tell, so comment out
  #CMC_rnaseq <- merge(CMC_rnaseq, CMC_rnaseq_synapse_table,by='Individual_ID', all=T)
  
  
  if(nrow(CMC_snp)!=length(unique(CMC_snp$Individual_ID)) | nrow(CMC_snp)!=length(unique(CMC_snp$Genotyping_Sample_ID))){ stop("More than 1 DNA per individual, change code")}
  CMC_snp[,c('260/280', '260/230') := NULL]
  CMC <- merge(CMC_rnaseq, CMC_snp, by=c('Individual_ID'), all.x=T)
  
  CMC <- merge(CMC, CMC_clinical, by='Individual_ID',all.x=T)

  CMC <- merge(CMC, CMC_dna_metadata, by=c('Individual_ID','Genotyping_Sample_ID','Genotyping_Chip_Type'),all.x=T)
  setnames(CMC, old=c('Genotyping_Sample_ID', 'Institution','Reported_Gender', 'Age_of_Death'), 
                new=c('genotype_id', 'cohort','Gender', 'age_death'))
  same_rows(CMC, list(CMC_rnaseq))
  
  CMC$SampleFull <- paste0(CMC$Individual_ID,'_',CMC$rnaseq_id)
  CMC$MetaCohort <- 'CMC'
  
  # remove some more columns that werent' removed before
  # better to only select the columns want to keep but like to be explicity in removing
  CMC[,c('SentrixBarcode','Genotyping_SentrixBarcode','Genotyping_SentrixPosition','DNA_isolation_Dneasy_Kit_ID','Date_Dissected',
               'DNA_isolation_Prep_Date','SCZ_Pair','BP_Pair','Flowcell_Batch','Library_Batch','Ribozero_Batch','RNA_isolation_Batch',
                'Study.x','Study.y','SentrixPosition','Dneasy_Kit_ID','DNA_Prep_Date','DNA_Prep_Operator', 'Dissection_Operator',
         'Flowcell','RNA_isolation_PrepDate','RNA_isolation_PrepOperator','DNA_isolation_Prep_Operator') := NULL]
  CMC$Brodmann_Area <- paste0('ba',CMC$Brodmann_Area)
  CMC[Brodmann_Area=='baNA', Brodmann_Area := NA]
  
  CMC[, Cause_of_Death:=as.character(Cause_of_Death)]
  CMC[Cause_of_Death=='1',Cause_of_Death:='Cardiovascular']
  CMC[Cause_of_Death=='2',Cause_of_Death:='Non-brain cancer']
  CMC[Cause_of_Death=='3',Cause_of_Death:='Inection and parasitic disease']
  CMC[Cause_of_Death=='4',Cause_of_Death:='COPD']
  CMC[Cause_of_Death=='5',Cause_of_Death:='Other']
  CMC[Cause_of_Death=='-1',Cause_of_Death:=NA]

  return(CMC)
}

parse_GTEx <- function(opt){
  GTEx <- fread(paste0(opt$metadataDir,'/GTEx/GTEx-RNASeq-samplesheet.txt'))
  GTEx <- GTEx[which(grepl('Brain',GTEx$`Comment[original body site annotation]`)),]
  GTEx <- GTEx[!duplicated(GTEx$`Comment[ENA_RUN]`),]
  
  setnames(GTEx, old=c('Comment[ENA_RUN]','Characteristics[individual]','Characteristics[sex]',
                       'Comment[INSTRUMENT_MODEL]','Comment[LIBRARY_LAYOUT]','Comment[LIBRARY_SELECTION]',
                       'FactorValue [organism part]','Comment[histological type]',
                       'Comment[original body site annotation]'),
                 new=c('rnaseq_id','genotype_id','Gender',
                       'instrument_model','library_layout','library_selection',
                       'organism_part','histological_part',
                       'body_site_annotations'))
  GTEx$MetaCohort <- 'GTEx'
  GTEx$cohort <- 'GTEx'
  # protocol REF is in twice, so remove already 1
  GTEx$`Protocol REF` <- NULL
  GTEx$`Protocol REF` <- NULL
  GTEx[,c('Comment[assembly]','Comment[SUBMITTED_FILE_NAME]','Scan Name','Comment[ENA_EXPERIMENT]','Technology type',
          'Comment[technical replicate group]','Comment[original technical replicate annotation]','Assay Name',
          'Performer','Comment[NOMINAL_LENGTH]','Comment[LIBRARY_STRATEGY]','Comment[LIBRARY_SOURCE]',
          'Extract Name','Protocol REF','Comment[gap_subject_id]','Comment[dbGAP_SAMPLE]','Comment[ENA_SAMPLE]',
          'Comment[BioSD_Sample]','Characteristics[clinical information]','Characteristics[disease]',
          'Characteristics[organism part]','Characteristics[organism]','histological_part'):=NULL]
  
  GTEx_sra_extra_samples <- fread(paste0(opt$metadataDir,'/GTEx/GTEx_SRA_samplesheet.txt'))
  GTEx <- rbind.fill(GTEx,GTEx_sra_extra_samples )
  GTEx$SampleFull <- paste0(GTEx$rnaseq_id,'_',GTEx$rnaseq_id)
  
  return(GTEx)
}

parse_NABEC <- function(opt){
  NABEC <- fread(paste0(opt$metadataDir,'/NABEC/NABEC_phenotypes.txt'))
  NABEC$individualID <- gsub('_e$|-wes\\.2$|-wes\\.1$|-cage-fctx$|-rna-fctx\\.1$|-rna-fctx\\.2$|-rna-fctx\\.3$|-rna-fctx\\.4$|-mrna-fctx$|-wes.3$|-wes.4$|-wes.5$',
                             '',NABEC$Sample_Name)
  
  NABEC_RNA <- NABEC[NABEC$Assay_Type=='RNA-Seq',]

  NABEC_WXS <- NABEC[NABEC$Assay_Type=='WXS',]
  
  NABEC_WXS$genotype_id <- NABEC_WXS$submitted_subject_id
  NABEC_WXS <- NABEC_WXS[,c('individualID','genotype_id')]
  
  

  NABEC <- merge(NABEC_RNA, NABEC_WXS, by='individualID', all.x=T)
  setnames(NABEC, old=c('sex','Center_Name','Run','LibraryLayout','LibrarySelection'), new=c('Gender','cohort','rnaseq_id','library_layout','library_selection'))
  # many columns are not informative, remove
  NABEC$SampleFull <- paste0(NABEC$BioSample,'_',NABEC$rnaseq_id)
  NABEC[,c('DATASTORE_filetype','DATASTORE_provider','Organism','Consent','LibrarySource','MBytes','study_design',
           'gap_parent_phs','biospecimen_repository','study_name','gap_accession','biospecimen_repository_sample_id',
           'analyte_type','Sample_Name','SRA_Study','SRA_Sample','MBases','Library_Name','Experiment','BioProject',
           'BioSample','Assay_Type','ReleaseDate','LoadDate'):=NULL]
  NABEC$MetaCohort <- 'NABEC'
  return(NABEC)
}

parse_TargetALS <- function(opt){
  #TODO: add more TargetALS samplesheets (check if more are available)
  TargetALS_metadata_rna <- fread(paste0(opt$metadataDir,'/TargetALS/Target_ALS_RNA_Metadata_03062019.txt'))
  TargetALS_metadata_rna$ExternalSampleId <- gsub('-2$','.2',TargetALS_metadata_rna$ExternalSampleId)
  TargetALS_metadata_rna$ExternalSampleId <- gsub('-','_',TargetALS_metadata_rna$ExternalSampleId)
  
  # to merge with wgs metadata, change to character
  setnames(TargetALS_metadata_rna, old='Age at Death', new='Age of Death')
  TargetALS_metadata_rna[, `Age of Death` := as.character(`Age of Death`)]
  
  TargetALS_metadata_wgs <- fread(paste0(opt$metadataDir,'/TargetALS/TALS_Data_Release_WGS_Metadata_10Oct2018.txt'))
  
  TargetALS_metadata_wgs$`Data File ID` <- NULL
  TargetALS_metadata_wgs <- fread(paste0(opt$metadataDir,'/TargetALS/TALS_Data_Release_WGS_Metadata_10Oct2018.txt'))
  TargetALS_metadata_wgs[TargetALS_metadata_wgs=='']<-NA
  TargetALS_metadata_wgs[`Age of Death`=='Unknown', `Age of Death` := NA]
  TargetALS_metadata_wgs[`Age of Death`=='90 or Older', `Age of Death` := '90+']
  TargetALS_metadata_wgs$`Sample Source` <- NULL
  TargetALS <- merge(TargetALS_metadata_rna,TargetALS_metadata_wgs,
                     by='ExternalSubjectId',all.x=T)
  
  # Remove columns with only NA
  TargetALS <- TargetALS[,which(unlist(lapply(TargetALS, function(x)!all(is.na(x))))),with=F]
  
  # Since some columns are in both RNAseq and WGS metadata, remove double columns and rename to remove the .x postfix
  drop.cols <- grep("\\.y$", colnames(TargetALS))
  TargetALS[, (drop.cols) := NULL]
  colnames(TargetALS) <- gsub('\\.x$','',colnames(TargetALS))
  

  # TargetALS RNAseq clincal file seems to only contain data already in the RNA metadata samplesheet. Will keep this
  # commented line in incase it's needed in the future
  # TargetALS_clinical1 <- fread(paste0(opt$metadataDir,'/TargetALS/UMG-Target_ALS_RNA_Clinical_Data_06122018.txt'))
  TargetALS_clinical2 <- fread(paste0(opt$metadataDir,'/TargetALS/UMG-Target_ALS_WGS_Clinical_Data_06122018.txt'))
  
  # Some values in TargetALS_clinical2 are empty, including for teh subject ID. Change these to NA and then remove all rows without subject ID
  TargetALS_clinical2[TargetALS_clinical2=='']<-NA
  TargetALS_clinical2 <- TargetALS_clinical2[!is.na(TargetALS_clinical2$ExternalSubjectId),]
  # remove columns with only NA
  TargetALS <- merge(TargetALS,TargetALS_clinical2, by='ExternalSubjectId',all.x=T)
  TargetALS <- TargetALS[!duplicated(TargetALS$ExternalSampleId),]
  TargetALS <- TargetALS[,which(unlist(lapply(TargetALS, function(x)!all(is.na(x))))),with=F]
  
  # Since some columns are in both TargetALS and clinical, remove double columns and rename to remove the .x postfix
  drop.cols <- grep("\\.y$", colnames(TargetALS))
  TargetALS[, (drop.cols) := NULL]
  colnames(TargetALS) <- gsub('\\.x$','',colnames(TargetALS))
  
  
  setnames(TargetALS, old=c('ExternalSampleId','CGND ID',       'Sex',  'Age of Death'), 
                    new=c('rnaseq_id',        'genotype_id','Gender','age_death'))
  TargetALS$SampleFull <- paste0(TargetALS$ExternalSubjectId,'_',TargetALS$rnaseq_id)
  TargetALS[,c('Quote', 'NeuroBank ID', 'Project','ExternalSubjectId',
               'Data File ID','NB Data downloaded') :=NULL]
  # some IDs have spaces in them, should be _
  TargetALS$MetaCohort <- 'TargetALS'
  #same_rows(TargetALS, list(TargetALS_metadata_rna))
  TargetALS$cohort <- 'TargetALS'
  TargetALS$rnaseq_id <- gsub('HRA-','HRA_',TargetALS$rnaseq_id)
  TargetALS$rnaseq_id <- gsub('CGND_','',TargetALS$rnaseq_id)
  TargetALS$rnaseq_id <- gsub('-2','',TargetALS$rnaseq_id)
  TargetALS$rnaseq_id <- gsub('\\.2','',TargetALS$rnaseq_id)
  same_rows(TargetALS_metadata_rna, list(TargetALS))
  return(TargetALS)
}

parse_ENA <- function(opt){
  ENA <- fread(paste0(opt$metadataDir,'/ENA/samplesheet_ENA_20181212.txt'))

  ENA <- ENA[,c('study_accession','sample_accession', 'run_accession', 'instrument_platform', 
                'instrument_model','library_layout','library_selection')]
  
  setnames(ENA, old=c('study_accession','sample_accession','run_accession'),new=c('cohort','genotype_id','rnaseq_id'))
  ENA$SampleFull <- paste0(ENA$genotype_id,'_',ENA$rnaseq_id)
  ENA$MetaCohort <- 'ENA'
  return(ENA)
}

parse_BrainGVEx <- function(opt){
  BrainGVEx_rnaseq_synapse_table <- fread(paste0(opt$metadataDir,'/BrainGVEx/BrainGVEx_RNA_table.txt'))
  BrainGVEx_rnaseq_synapse_table <-   BrainGVEx_rnaseq_synapse_table[,which(unlist(lapply(BrainGVEx_rnaseq_synapse_table, function(x)!all(is.na(x))))),with=F]
  BrainGVEx_rnaseq_synapse_table[,c('id','name','consortium','species','study','parentId','Contributor','grant',
                           'dataFileHandleId','currentVersion','createdBy','Capstone_1','assay','organ',
                           'Capstone_4','nucleicAcidSource','fileFormat','isStranded','readLength','runType','BrodmannArea') := NULL]
  BrainGVEx_rnaseq_metadata <- fread(paste0(opt$metadataDir,'/BrainGVEx/UIC-UChicago-U01MH103340_BrainGVEX_RNAseqMetadata_August2016Release.csv'))
  BrainGVEx_rnaseq_metadata[,c('rRNARate','MappedReads_Multimapped','MappedReads_Primary','TotalReads','ReadLength','LibraryKit','LibraryBatch','RNAIsolationBatch',
                             'CellType','File_Name','Assay') := NULL]
  setnames(BrainGVEx_rnaseq_metadata, old=c('Individual_ID (RNAseq Library BID, Synapse)'), new =c('individualID'))
  BrainGVEx_rnaseq_metadata <- BrainGVEx_rnaseq_metadata[BrainGVEx_rnaseq_metadata$individualID != '',]
  
  BrainGVEx_rnaseq <- merge(BrainGVEx_rnaseq_synapse_table, BrainGVEx_rnaseq_metadata, by='individualID',all.x=T)
  
  BrainGVEx_clinical <- fread(paste0(opt$metadataDir,'/BrainGVEx/UIC-UChicago-U01MH103340_BrainGVEX_ClinicalMetadata_AgeCensored_August2016Release.tsv'))

  BrainGVEx_clinical <-   BrainGVEx_clinical[,which(unlist(lapply(BrainGVEx_clinical, function(x)!all(is.na(x))))),with=F]
  
  BrainGVEx_clinical[,c('GenotypeData','Opiates','Cocaine','THC','Nic_Cot','MoodStab','Benzos','Organism','Grant','index','StudyName') := NULL]
  setnames(BrainGVEx_clinical, old=c('Individual_ID (RNAseq library BID, Primary ID)'), new =c('individualID'))
  
  BrainGVEx <- merge(BrainGVEx_rnaseq, BrainGVEx_clinical, by='individualID',all.x=T)
  setnames(BrainGVEx, old=c('specimenID'), new=c('rnaseq_id'))
  BrainGVEx$cohort <- 'BrainGVEx'
  BrainGVEx$MetaCohort <- 'PsychEncode'
  return(BrainGVEx)
}

parse_UCLA_ASD <- function(opt){
  rna_metadata <- fread(paste0(opt$metadataDir,'/UCLA_ASD/UCLA_R01MH094714_ASD_Metadata_RNAseq.csv'))
  rna_metadata <- rna_metadata[,c('Individual_ID','Sample_ID','BrodmannArea','BrainRegion','RIN','LibraryPrep','LibraryKit','RunType','SequencingPlatform')]
  snp_metadata <- fread(paste0(opt$metadataDir,'/UCLA_ASD/KCL_R01MH094714_ASD_Metadata_Illumina450K.csv'))
  snp_metadata <- snp_metadata[,c('Individual_ID','BroadmannArea','BrainRegion','Diagnosis','Sex','Age','PMI','BrainBank','PrimaryCauseOfDeath','CETs','DNAmAge')]
  clinical_metadata <- fread(paste0(opt$metadataDir,'/UCLA_ASD/UCLA_R01MH094714_ASD_Metadata_Clinical_August2016Release.csv'))
  clinical_metadata <- clinical_metadata[,c('Grant','StudyName','Organism'):=NULL]
  
  
  UCLA_ASD <- merge(rna_metadata, snp_metadata, by=c('Individual_ID','BrainRegion'),all.x=T,allow.cartesian=TRUE)
  UCLA_ASD <- merge(UCLA_ASD, clinical_metadata, by='Individual_ID',all.x=T,allow.cartesian=TRUE)
    
  setnames(UCLA_ASD, old=c('Individual_ID','Sample_ID'), new = c('individualID','rnaseq_id'))
  
  UCLA_ASD[,rnaseq_id := gsub('Sample_','',rnaseq_id)]
  UCLA_ASD[,rnaseq_id := gsub('_3rd','',rnaseq_id)]
  UCLA_ASD[,rnaseq_id := gsub('_2nd','',rnaseq_id)]
  UCLA_ASD[,rnaseq_id := gsub('_1st','',rnaseq_id)]
  
  UCLA_ASD[rnaseq_id=='UMB5176_ba41-42-22',rnaseq_id:= 'UMB5176_ba41.42.22'] 
  UCLA_ASD[rnaseq_id=='AN03345_ba41-42-22',rnaseq_id:= 'AN03345_ba41.42.22'] 
  UCLA_ASD[rnaseq_id=='UMB1376_ba41-42-22',rnaseq_id:= 'UMB1376_ba41.42.22'] 
  UCLA_ASD[rnaseq_id=='AN15088_ba41-42-22',rnaseq_id:= 'AN15088_ba41.42.22'] 
  UCLA_ASD[rnaseq_id=='UMB4337_ba41-42-22',rnaseq_id:= 'UMB4337_ba41.42.22'] 
  UCLA_ASD[rnaseq_id=='AN11864_ba41-42-22',rnaseq_id:= 'AN11864_ba41.42.22'] 
  UCLA_ASD[rnaseq_id==']AN04479_BA7',rnaseq_id:= 'AN04479_BA7'] 
  
  
  same_rows(rna_metadata, list(UCLA_ASD))
  UCLA_ASD$cohort <- 'UCLA_ASD'
  UCLA_ASD$MetaCohort <- 'PsychEncode'
  return(UCLA_ASD)
}

parse_BipSeq <- function(opt){
  BipSeq_rnaseq_metadata <- fread(paste0(opt$metadataDir,'/BipSeq/LIBD-JHU_BipSeq_R01MH105898_Metadata_RNAseq_DLPFC_Limbic.csv'))
  BipSeq_rnaseq_metadata <- BipSeq_rnaseq_metadata[,c('individualID','individualIdSource','specimenID','tissue','brodmannArea','cellType','RIN','platform','libraryPrep','runType','isStranded')]
  BipSeq_rnaseq_metadata <- BipSeq_rnaseq_metadata[!BipSeq_rnaseq_metadata$individualID == "",]
  BipSeq_snp_metadata <- fread(paste0(opt$metadataDir,'/BipSeq/LIBD-JHU_R01MH105898_BipSeq_Metadata_SNPgenotypes_August2016Release.csv'))
  BipSeq_snp_metadata <- BipSeq_snp_metadata[,c('Individual_ID','platform')]
  BipSeq_clinical_metadata <- fread(paste0(opt$metadataDir,'/BipSeq/LIBD-JHU_R01MH105898_BipSeq_Metadata_Clinical.csv'))
  BipSeq_clinical_metadata <- BipSeq_clinical_metadata[,c('grant','study','Organism'):=NULL]
  
  setnames(BipSeq_snp_metadata, old=('platform'),new='SnpArrayPlatform')
  BipSeq <- merge(BipSeq_rnaseq_metadata, BipSeq_snp_metadata, by.x='individualID',by.y='Individual_ID',all.x=T)
  BipSeq <- merge(BipSeq, BipSeq_clinical_metadata, by.x='individualID',by.y='Individual_ID',all.x=T)
  setnames(BipSeq, old=c('individualID','specimenID'), new = c('individualID','rnaseq_id'))
  BipSeq$cohort <- 'BipSeq'
  BipSeq$MetaCohort <- 'PsychEncode'
  
  return(BipSeq)
}

parse_CMC_HBCC <- function(opt){
  CMC_HBCC_rnaseq <- fread(paste0(opt$metadataDir,'/CMC_HBCC/CMC_HBCC_RNA_table.txt'))
  CMC_HBCC_rnaseq <- CMC_HBCC_rnaseq[,c('individualID','specimenID','organ','tissue','BrodmannArea','nucleicAcidSource','hemisphere','PMI','pH','libraryPrep',
                                     'RIN','isStranded','platform')]
  setnames(CMC_HBCC_rnaseq, old=c('specimenID'),new=c('rnaseq_id'))
  CMC_HBCC_rnaseq$genotype_id <- CMC_HBCC_rnaseq$individualID
  CMC_HBCC_rnaseq$cohort <- 'CMC_HBCC'
  CMC_HBCC_rnaseq$MetaCohort <- 'PsychEncode'
  
  return(CMC_HBCC_rnaseq)
}

