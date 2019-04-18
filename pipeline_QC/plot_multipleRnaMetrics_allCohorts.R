
library(data.table)
library(ggplot2)
library(jcolors)
library(plyr)
library("dplyr")
library(stringr)
library("optparse")
library(RColorBrewer)

# Get command line arguments 
option_list = list(
  make_option(c("-i", "--input"), type="character", default=getwd(), 
              help="path to input file base directory (will recursive search)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir", metavar="character"),
  make_option(c("-p", "--pcaOutliers"), type="character", default=NULL,
              help="path to file with PCA outliers", metavar="character"),
  make_option(c("-n", "--NABEC_pheno"), type="character", default=NULL,
              help="phenotype data for NABEC", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

test <- function(){
  opt <- list()
  opt$input <- "/Users/NPK/UMCG/projects/biogen/cohorts/"
  opt$output <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/QC"
  opt$pcaOutliers <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/QC/rna-pcaoutliers.txt"
  opt$NABEC_pheno <- "/Users/NPK/UMCG/projects/biogen/cohorts/NABEC/phenotypes/NABEC_phenotypes.txt"
  return(opt)
}
# comment out when not testing
opt <- test()


print(paste('Searching for input files in: ', opt$input))
if(!is.null(opt$pcaOutliers)){
  print(paste('Reading pca outlier file from: ',opt$pcaOutliers))
}
print(paste('Writing output files to: ',opt$output))


##### read in samples that are filtered out during PC #####
if(!is.null(opt$pcaOutliers)){
    pca_filtered_samples <- fread(opt$pcaOutliers,header=F)
}
#####
nabec_phenotypes <- NULL
#### Read in NABEC phenotype data ####
if(!is.null(opt$nabec_phenotypes)){
  nabec_phenotypes <- fread(opt$NABEC_pheno)
}
####

##### get cohort from path function that's being used a couple of times ####
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
  stop(paste('unknown cohort:',path))
}
#####


##### Read RnaMetrics alignment QC #######
print('Readin RNAseq QC metrics files')
RnaMetrics_multiQC_files <- list.files(path=opt$input, 
                                       pattern = "multiqc_picard_RnaSeqMetrics.txt$", 
                                       recursive = TRUE,
                                       full.names = T)
BPD <- RnaMetrics_multiQC_files[grepl("BPD", RnaMetrics_multiQC_files)]

RnaMetrics_multiQC_files <- RnaMetrics_multiQC_files[!grepl("BPD", RnaMetrics_multiQC_files)]
RnaMetric_qc_all <- data.frame()
for(f in RnaMetrics_multiQC_files){
  cohort <- get_cohort_from_path(f)
  if(cohort == 'ENA'){
    if(!grepl('multiQC_alignment',  f)) { next }
  }
  else{print(f)
    print(cohort)}
  
  RnaMetrics_qc <- fread(f)
  RnaMetrics_qc$cohort <- cohort
  RnaMetrics_qc$SAMPLE <- NULL
  RnaMetrics_qc$study <- NA
  if(cohort == 'ENA'){
    pattern <- ".*/data/(.*)_batch[0-9]+.*"
    RnaMetrics_qc$study <- sub(pattern, "\\1", f)
  }
  RnaMetric_qc_all <- rbind(RnaMetric_qc_all, RnaMetrics_qc)
}
RnaMetric_qc_all <- data.frame(RnaMetric_qc_all)
RnaMetric_qc_all$Sample <- gsub('.cram','',RnaMetric_qc_all$Sample)
RnaMetric_qc_all$Sample <- gsub('-','_',RnaMetric_qc_all$Sample)

if('TargetALS' %in% RnaMetric_qc_all$cohort){
  RnaMetric_qc_all[RnaMetric_qc_all$cohort=="TargetALS",]$Sample <- str_match(RnaMetric_qc_all[RnaMetric_qc_all$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
}
######

###### Read MultiMetrics QC #####
print('Readin MultipleMetrics files')
MultipleMetrics_multiQC_files <- list.files(path=opt$input, 
                                            pattern = "multiqc_picard_AlignmentSummaryMetrics.txt$", 
                                            recursive = TRUE,
                                            full.names = T)

MultipleMetrics_multiQC_files <- MultipleMetrics_multiQC_files[!grepl("BPD", MultipleMetrics_multiQC_files)]

MultipleMetric_qc_all <- data.frame()
for(f in MultipleMetrics_multiQC_files){
  cohort <- get_cohort_from_path(f)
  if(cohort == 'ENA'){
    if(!grepl('multiQC_alignment',  f)) { next }
  }
  MultipleMetrics_qc <- fread(f)
  MultipleMetrics_qc$cohort <- cohort
  MultipleMetrics_qc$SAMPLE <- NULL
  MultipleMetric_qc_all <- rbind(MultipleMetric_qc_all, MultipleMetrics_qc)
}
MultipleMetric_qc_all <- data.frame(MultipleMetric_qc_all)
MultipleMetric_qc_all$Sample <- gsub('.cram','',MultipleMetric_qc_all$Sample)
MultipleMetric_qc_all$Sample <- gsub('-','_',MultipleMetric_qc_all$Sample)

if('TargetALS' %in% MultipleMetric_qc_all$cohort){
  MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="TargetALS",]$Sample <- str_match(MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
}
####### 

###### Read star QC #####
print('Reading star QC files')
STAR_multiQC_files <- list.files(path=opt$input, 
                                 pattern = "multiqc_star.txt$", 
                                 recursive = TRUE,
                                 full.names = T)
STAR_multiQC_files <- STAR_multiQC_files[!grepl("BPD", STAR_multiQC_files)]

STAR_qc_all <- data.frame()
for(f in STAR_multiQC_files){
  STAR_qc <- fread(f)
  STAR_qc$cohort <- get_cohort_from_path(f)
  STAR_qc_all <- rbind(STAR_qc_all, STAR_qc)
}
STAR_qc_all <- data.frame(STAR_qc_all)
STAR_qc_all$Sample <- gsub('.cram','',STAR_qc_all$Sample)
STAR_qc_all$Sample <- gsub('-','_',STAR_qc_all$Sample)

if('TargetALS' %in% STAR_qc_all$cohort){
  STAR_qc_all[STAR_qc_all$cohort=="TargetALS",]$Sample <- str_match(STAR_qc_all[STAR_qc_all$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
}
####### 

###### Read fastqc QC #####
print('Readin fastqc files')
FastQC_multiQC_files <- list.files(path=opt$input, 
                                   pattern = "multiqc_fastqc.txt$", 
                                   recursive = TRUE,
                                   full.names=T)
FastQC_multiQC_files <- FastQC_multiQC_files[!grepl("BPD", FastQC_multiQC_files)]

FastQC_qc_all <- data.frame()
for(f in FastQC_multiQC_files){
  FastQC_qc <- fread(f)

  FastQC_qc$cohort <- get_cohort_from_path(f)
  FastQC_qc_all <- rbind(FastQC_qc_all, FastQC_qc, fill=T)
}
# Because FastQC gives results per fastq file instead of per sample, adjust the names and merge them together
FastQC_qc_all <- data.frame(FastQC_qc_all)
FastQC_qc_all$Sample <- gsub('.cram','',FastQC_qc_all$Sample)

FastQC_qc_all_R1 <- FastQC_qc_all[grepl('*_R1$|*\\.R1$|*_1$', FastQC_qc_all$Sample),]
FastQC_qc_all_R1$Sample <- gsub('*_R1$|*\\.R1$|*_1$','',FastQC_qc_all_R1$Sample)
FastQC_qc_all_R2 <- FastQC_qc_all[grepl('*_R2$|*\\.R2$|*_2$', FastQC_qc_all$Sample),]
FastQC_qc_all_R2$Sample <- gsub('*_R2$|*\\.R2$|*_2$','',FastQC_qc_all_R2$Sample)

# everything that is leftover is single end, save as R1
single_end <- FastQC_qc_all[!grepl('*_R1$|*\\.R1$|*_1$', FastQC_qc_all$Sample),]
single_end <- single_end[!grepl('*_R2$|*\\.R2$|*_2$', single_end$Sample),]
FastQC_qc_all_R1 <- rbind(FastQC_qc_all_R1, single_end)

# fill.x = True because single end sample are only in R1
FastQC_qc_all_merged <- merge(FastQC_qc_all_R1, FastQC_qc_all_R2, by=c('Sample','cohort'), suffix = c("_R1", "_R2"),
                              all.x=T)
FastQC_qc_all_merged$FastQC_original_sample <- FastQC_qc_all_merged$Sample

# For a few samples fastqc also got information on a 3rd file that contained orphan reads. Remove these (recognizable by lower total reads
FastQC_qc_all_merged <- FastQC_qc_all_merged[order(FastQC_qc_all_merged$Sample, -abs(FastQC_qc_all_merged$Total.Sequences_R1) ), ] #sort by id and reverse of abs(value)
FastQC_qc_all_merged <- FastQC_qc_all_merged[ !duplicated(FastQC_qc_all_merged$Sample), ]    

# TargetALS sample naming is different between FastQC and picard results, change them here to make them match later
if('TargetALS' %in% FastQC_qc_all_merged$cohort){
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- gsub('-','_',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample )
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- str_match(FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample , ".*(HRA_[0-9]+)")[, 2]
    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- gsub('-b38','',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample)
}
# Same for Brainseq
if('Brainseq' %in% FastQC_qc_all_merged$cohort){
    Brainseq_full_sample <- MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="Brainseq",]$Sample
    Brainseq_fastqc_sample <- FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="Brainseq",]$Sample

    FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="Brainseq",]$Sample <- unlist(
                                                                            lapply(Brainseq_fastqc_sample, function(x) {
                                                                                Brainseq_full_sample[grepl(x, Brainseq_full_sample)]
                                                                            }))
}


####### 



###### merge FastQC and picard metrics together #####
RnaMetric_qc_all[grepl('PCT', colnames(RnaMetric_qc_all))] <- RnaMetric_qc_all[grepl('PCT', 
                                                                                colnames(RnaMetric_qc_all))]/100






CombinedMetrics <- merge(RnaMetric_qc_all, MultipleMetric_qc_all, by=c('Sample','cohort',
                                                                       'PF_ALIGNED_BASES',
                                                                       'LIBRARY',
                                                                       'READ_GROUP','PF_ALIGNED_BASES'),fill=T)
CombinedMetrics <- merge(CombinedMetrics, STAR_qc_all, by=c('Sample','cohort'), all = TRUE)
CombinedMetrics$SampleFull <- CombinedMetrics$Sample

if('Braineac' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="Braineac",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="Braineac",]$Sample, ".*(A653.*)")[, 2]
}
if('GTEx' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="GTEx",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="GTEx",]$Sample, "(.*)_.*")[, 2]
}
if('NABEC' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="NABEC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="NABEC",]$Sample, ".*_(.*)")[, 2]
}
if('ENA' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="ENA",]$Sample <- str_match(
                                                                CombinedMetrics[CombinedMetrics$cohort=="ENA",]$Sample,
                                                                                                         ".*_(.*)")[, 2]
}
if('TargetALS' %in% CombinedMetrics$cohort){
    CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample, ".*(HRA.*)")[, 2]
    CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample <- gsub('-2','',CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample)
}

if('BrainGVEx' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample <- gsub('_','-',CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample )
  CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="BrainGVEx",]$Sample,"^([0-9]+-[0-9]+)-[0-9]+-[0-9]+")[, 2]
}
if('BipSeq' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="BipSeq",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="BipSeq",]$Sample,"Br[0-9]+_(R[0-9]+)")[, 2]
}
if('UCLA_ASD' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="UCLA_ASD",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="UCLA_ASD",]$Sample,"[aA-zZ]+[0-9]+_(.+)")[, 2]
  
}
if('CMC_HBCC' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="CMC_HBCC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="CMC_HBCC",]$Sample,"individualID.*_specimenID.(.*)")[, 2]
  
}
if('CMC' %in% CombinedMetrics$cohort){
  CombinedMetrics[CombinedMetrics$cohort=="CMC",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="CMC",]$Sample,"individualID.*_specimenID.(.*)")[, 2]
}


CombinedMetricsWithFastQC <- merge(CombinedMetrics, FastQC_qc_all_merged, by=c('Sample','cohort'), all = TRUE)
######


##### read in AMP-D ######
AMP_AD_multiMetrics_files <- list.files(path=opt$input,
                                        pattern = "CombinedMetrics.csv$", 
                                        recursive = T,
                                        full.names = T)
if(length(AMP_AD_multiMetrics_files) > 0){
    AMP_AD_multimetrics_qc_all <- data.frame()
    for(f in AMP_AD_multiMetrics_files){
      multimetrics_qc <- fread(f)
      multimetrics_qc$cohort <- sapply(strsplit(basename(f), "_CombinedMetrics"), "[[", 1)
      AMP_AD_multimetrics_qc_all <- rbind(AMP_AD_multimetrics_qc_all, multimetrics_qc)
    }
    AMP_AD_multimetrics_qc_all <- data.frame(AMP_AD_multimetrics_qc_all)
    colnames(AMP_AD_multimetrics_qc_all) <- gsub('.*__','',colnames(AMP_AD_multimetrics_qc_all))
    names(AMP_AD_multimetrics_qc_all)[names(AMP_AD_multimetrics_qc_all) == 'sample'] <- 'Sample'
    ######

    ####### merge AMP-AD with the rest #######
    AMP_AD_multimetrics_qc_all$PF_ALIGNED_BASES.1 <- NULL
    all_cohorts <- rbind.fill(CombinedMetricsWithFastQC, AMP_AD_multimetrics_qc_all)
}else{
    all_cohorts <- CombinedMetricsWithFastQC
}
all_cohorts <- all_cohorts[order(all_cohorts$cohort),]
all_cohorts$sampleID <- 1:nrow(all_cohorts)
#######

###### columns with % dont actually show %, have to multiply with 100
all_cohorts[which(grepl('PCT',colnames(all_cohorts)))] <- all_cohorts[which(grepl('PCT',colnames(all_cohorts)))]*100

###### set filter based on 3 measurements ######
all_cohorts$FILTER <- "NO"
print('set filter')


if(sum(!is.na(all_cohorts$PCT_CODING_BASES) & all_cohorts$PCT_CODING_BASES <= 10) > 0){
    all_cohorts[!is.na(all_cohorts$PCT_CODING_BASES) & all_cohorts$PCT_CODING_BASES <= 10,]$FILTER <- "YES"
}
if(sum(!is.na(all_cohorts$PCT_PF_READS_ALIGNED) & all_cohorts$PCT_PF_READS_ALIGNED <= 60) > 0){
    all_cohorts[!is.na(all_cohorts$PCT_PF_READS_ALIGNED) & all_cohorts$PCT_PF_READS_ALIGNED <= 60,]$FILTER <- "YES"
}


# these are wrong in the expression table, but since we don't want to rewrite that one change it here
all_cohorts[all_cohorts$Sample=='UMB5176_ba41_42_22',]$Sample <- 'UMB5176_ba41.42.22'
all_cohorts[all_cohorts$Sample=='AN03345_ba41_42_22',]$Sample <- 'AN03345_ba41.42.22'
all_cohorts[all_cohorts$Sample=='UMB1376_ba41_42_22',]$Sample <- 'UMB1376_ba41.42.22'
all_cohorts[all_cohorts$Sample=='AN15088_ba41_42_22',]$Sample <- 'AN15088_ba41.42.22'
all_cohorts[all_cohorts$Sample=='UMB4337_ba41_42_22',]$Sample <- 'UMB4337_ba41.42.22'
all_cohorts[all_cohorts$Sample=='AN11864_ba41_42_22',]$Sample <- 'AN11864_ba41.42.22'

# sample should have been merged but resequenced sample didn't get transfered to us until the first one was already processed.
# only keep cov. data for teh correct one

all_cohorts <- all_cohorts[!duplicated(all_cohorts$Sample),]



all_cohorts[all_cohorts$Sample=='AN04479_BA7',]


filtered_samples_file <- paste0(opt$output,"/samples_to_filter.txt")
write.table(all_cohorts[all_cohorts$FILTER=='YES',]$Sample, filtered_samples_file, quote=F, sep="\t", row.names=F)
print(paste0("Written samples to filer to ",filtered_samples_file))


nSamples <- all_cohorts %>% group_by(cohort) %>% dplyr::summarize(no = sum(FILTER=="NO"),yes = sum(FILTER=="YES"))
nSamples$total <- nSamples$no + nSamples$yes

all_cohorts$SampleFull <- gsub('individualID.','', all_cohorts$SampleFull)
all_cohorts$SampleFull <- gsub('specimenID.','', all_cohorts$SampleFull)

all_cohorts$pcaFiltered <- 'NO'

if(!is.null(opt$pcaOutliers)){
    # Add a separate column for those that got filtered with PCA
    all_cohorts[all_cohorts$SampleFull %in% pca_filtered_samples$V1,]$pcaFiltered <- 'YES'
}
#####

##### Merge NABEC phenotype into all_cohorts
all_cohorts$lib_selection <- 'unknown'
if(!is.null(nabec_phenotypes)){
  if('NABEC' %in% all_cohorts$cohort){
    nabec_phenotypes$Sample <- paste0(nabec_phenotypes$BioSample, '_', nabec_phenotypes$Run)
    all_cohorts[all_cohorts$cohort=="NABEC",]$lib_selection <- nabec_phenotypes[match(all_cohorts[all_cohorts$cohort=="NABEC",]$Sample,
                                                                             nabec_phenotypes$Sample),]$LibrarySelection
  }
}

#####

##### plot function so that it can be more easily repeated
QC_plotter <- function(all_cohorts, column, plot_pca_outliers = F, single_cohort =F, lib_selection=F){
  # add the super small nuber in case y-axis needs to be log scaled, but only for those columns where the max value > 100000 (so that columns with e.g. percentages don't get log scaled)
  if(max(all_cohorts[!is.na(all_cohorts[column]),][column])> 10000){
    all_cohorts[!is.na(all_cohorts[column]),][column] <-all_cohorts[!is.na(all_cohorts[column]),][column] +0.00000000000000000000001
  }
  
  # Because we have more than 8 colours, have to extend the Dark2 colour palette
  colourCount = length(unique(all_cohorts$cohort))
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  dark2ExtendedPallete <- getPalette(colourCount)
  
  p <- ggplot()
  
  if(single_cohort){
    p <- p + geom_point(data=all_cohorts, aes_string('sampleID', column, colour='FILTER'))
  }else{
    # if PCA outliers are plotted, make all others have lower alpha settings
    if(plot_pca_outliers){
      if(lib_selection){
        p <- p + geom_point(alpha=0.1)+ 
          geom_point(data=all_cohorts[all_cohorts$pcaFiltered=='YES',],  aes_string('sampleID', column, colour='cohort', shape='lib_selection'))
      }else{
        p <- p + geom_point(alpha=0.1)+ 
             geom_point(data=all_cohorts[all_cohorts$pcaFiltered=='YES',],  aes_string('sampleID', column, colour='cohort', shape='FILTER'))
      }
    }else{
      if(lib_selection){
        p <- p + geom_point(data=all_cohorts, aes_string('sampleID', column, colour='cohort', shape='lib_selection'))
      }else{
        p <- p + geom_point(data=all_cohorts, aes_string('sampleID', column, colour='cohort', shape='FILTER'))
      }
    }
  }
  p <- p + theme_bw(base_size=18)+
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())+  
    xlab('Samples')#+
    #facet_wrap(~cohort, scale='free_x')
  
  if(single_cohort){
    p <- p+scale_colour_brewer(palette="Dark2")
  }else{
    p <- p+scale_colour_manual(values=dark2ExtendedPallete)
  }

  # Add the log scaling for columns with high values
  if(max(all_cohorts[!is.na(all_cohorts[column]),][column])> 10000){
    p <- p + scale_y_continuous(trans='log10') +
      ylab(paste0('log10( ',column, ' )'))
  }
  
  # Plot a line where we put the threshold
  if(column == "PCT_CODING_BASES"){
    p <- p + geom_hline(yintercept=10, colour="red",  linetype="dashed")
  }
  if(column == "PCT_PF_READS_ALIGNED"){
    p <- p + geom_hline(yintercept=60, colour="red",  linetype="dashed")
  }
  return(p)
}
#####

##### plot STAR + multimetrics together #####

# make sure that the columns used for filtering are plotted first so that they are easy to find
columns_to_plot <- colnames(select_if(all_cohorts, is.numeric))
filter_columns <-  c("PCT_CODING_BASES","PCT_PF_READS_ALIGNED")
columns_to_plot <- columns_to_plot[!columns_to_plot %in% filter_columns]
columns_to_plot <- c(filter_columns, columns_to_plot)

print('start plotting')
pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots.pdf'), width=8, height=8)
for(column in columns_to_plot){
  if(column == "sampleID"){
    next
  }
  print(column)
  p <- QC_plotter(all_cohorts, column, FALSE)
  ggsave(paste0(opt$output,'/figures/QC_figures_separate/',column,'.png'), width=12, height = 8, plot=p)
  print(p)
}
dev.off()

colourCount = length(unique(all_cohorts$cohort))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
dark2ExtendedPallete <- getPalette(colourCount)
ggplot(all_cohorts, aes(cohort,fill=cohort,colour='yes'))+
  geom_bar()+ 
  theme_bw(base_size=18)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F)+ 
  coord_flip()+
  guides(fill=guide_legend(ncol=2))
ggsave(paste0(opt$output,'/figures/n_RNAseq_samples.png'),width=12, height=8)
#####

##### plot lib select #####
plot_with_lib_selection <- function(){
  # used for debugging nabec, not necesarry to run
  print('plot library selection')
  pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots_library_selection.pdf'), width=8, height=8)
  for(column in columns_to_plot){
    if(column == "sampleID"){
      next
    }
    print(column)
    p <- QC_plotter(all_cohorts, column, FALSE, lib_selection=T)
    print(p)
  }
  dev.off()
}
# plot_with_lib_selection()
#####

##### plot STAR + multimetrics again, but highlight PCA filtered samples #####
plot_with_pca_outliers <- function(){
  if(!is.null(opt$pcaOutliers)){
      print('plotting with pca outliers')
      pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots.highlightPcaFilteredSamples.pdf'), width=8, height=8)
      for(column in columns_to_plot){
        if(column == "sampleID"){
          next
        }
        print(column)
        p <- QC_plotter(all_cohorts, column, TRUE)
        ggsave(paste0(opt$output,'/figures/QC_figures_separate/',column,'.highlightPcaFilteredSamples.png'), width=12, height = 8, plot=p)
        print(p)
      }
      dev.off()
  }
}
#plot_with_pca_outliers()
#####
write.table(all_cohorts,file=paste0(opt$output,'/2019-04-09-Freeze2.TMM.Covariates.txt'),quote=F, row.names=T, sep='\t')
