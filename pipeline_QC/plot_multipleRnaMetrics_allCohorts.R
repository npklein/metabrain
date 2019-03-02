
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
              help="path to file with PCA outliers", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(paste('Searching for input files in: ', opt$input))
if(!is.null(opt$pcaOutliers)){
  print(paste('Reading pca outlier file from: ',opt$pcaOutliers))
}
print(paste('Writing output files to: ',opt$output))

#opt$input <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/QC"
#opt$output <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/QC"
#opt$pcaOutliers <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/QC/rna-pcaoutliers.txt"

##### read in samples that are filtered out during PC #####
if(!is.null(opt$pcaOutliers)){
    pca_filtered_samples <- fread(opt$pcaOutliers,header=F)
}
#####

##### get cohort from path function that's being used a couple of times ####
get_cohort_from_path <- function(path){
  # Below only works because all files are in same directory structure:
  # /path/to/<cohort_name>/QC/multiqc_data/file.txt
  # This extracts <cohort_name>
  cohort <- basename(dirname(dirname(dirname(f))))
  return(cohort)
}
#####


##### Read RnaMetrics alignment QC #######
RnaMetrics_multiQC_files <- list.files(path=opt$input, 
                                       pattern = "multiqc_picard_RnaSeqMetrics.txt$", 
                                       recursive = TRUE,
                                       full.names = T)
RnaMetrics_multiQC_files <- RnaMetrics_multiQC_files[!grepl("BPD", RnaMetrics_multiQC_files)]
RnaMetric_qc_all <- data.frame()
for(f in RnaMetrics_multiQC_files){
  RnaMetrics_qc <- fread(f)
  RnaMetrics_qc$cohort <- get_cohort_from_path(f)
  RnaMetrics_qc$SAMPLE <- NULL
  RnaMetric_qc_all <- rbind(RnaMetric_qc_all, RnaMetrics_qc)
}
RnaMetric_qc_all <- data.frame(RnaMetric_qc_all)
RnaMetric_qc_all$Sample <- gsub('.cram','',RnaMetric_qc_all$Sample)
RnaMetric_qc_all$Sample <- gsub('-','_',RnaMetric_qc_all$Sample)
######

###### Read MultiMetrics QC #####
MultipleMetrics_multiQC_files <- list.files(path=opt$input, 
                                            pattern = "multiqc_picard_AlignmentSummaryMetrics.txt$", 
                                            recursive = TRUE,
                                            full.names = T)
MultipleMetrics_multiQC_files <- MultipleMetrics_multiQC_files[!grepl("BPD", MultipleMetrics_multiQC_files)]

MultipleMetric_qc_all <- data.frame()
for(f in MultipleMetrics_multiQC_files){
  MultipleMetrics_qc <- fread(f)
  MultipleMetrics_qc$cohort <- get_cohort_from_path(f)
  MultipleMetrics_qc$SAMPLE <- NULL
  MultipleMetric_qc_all <- rbind(MultipleMetric_qc_all, MultipleMetrics_qc)
}
MultipleMetric_qc_all <- data.frame(MultipleMetric_qc_all)
MultipleMetric_qc_all$Sample <- gsub('.cram','',MultipleMetric_qc_all$Sample)
MultipleMetric_qc_all$Sample <- gsub('-','_',MultipleMetric_qc_all$Sample)
####### 

###### Read star QC #####
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
####### 

###### Read fastqc QC #####
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

FastQC_qc_all_R1 <- FastQC_qc_all[grepl('*\\.R1|*_1', FastQC_qc_all$Sample),]
FastQC_qc_all_R1$Sample <- gsub('.R1','',FastQC_qc_all_R1$Sample)
FastQC_qc_all_R1[FastQC_qc_all_R1$cohort=="GTEx",]$Sample <- gsub('_1','',FastQC_qc_all_R1[FastQC_qc_all_R1$cohort=="GTEx",]$Sample)
FastQC_qc_all_R2 <- FastQC_qc_all[grepl('*\\.R2|*_2', FastQC_qc_all$Sample),]
FastQC_qc_all_R2$Sample <- gsub('.R2','',FastQC_qc_all_R2$Sample)
FastQC_qc_all_R2[FastQC_qc_all_R2$cohort=="GTEx",]$Sample <- gsub('_2','',FastQC_qc_all_R2[FastQC_qc_all_R2$cohort=="GTEx",]$Sample)
FastQC_qc_all_merged <- merge(FastQC_qc_all_R1, FastQC_qc_all_R2, by=c('Sample','cohort'), suffix = c("_R1", "_R2"))
FastQC_qc_all_merged$FastQC_original_sample <- FastQC_qc_all_merged$Sample


# TargetALS sample naming is different between FastQC and picard results, change them here to make them match later
FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- str_match(FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample , ".*(HRA.*?)_.*")[, 2]
FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- gsub('-b38','',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample)
FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample <- gsub('-','_',FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="TargetALS",]$Sample )

# Same for Brainseq
Brainseq_full_sample <- MultipleMetric_qc_all[MultipleMetric_qc_all$cohort=="Brainseq",]$Sample
Brainseq_fastqc_sample <- FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="Brainseq",]$Sample

#FastQC_qc_all_merged[FastQC_qc_all_merged$cohort=="Brainseq",]$Sample  <- 
  
  t <- unlist(
                                                                            lapply(Brainseq_fastqc_sample, function(x) {
                                                                                Brainseq_full_sample[grepl(x, Brainseq_full_sample)]
                                                                            }))
####### 



###### merge FastQC and picard metrics together #####
RnaMetric_qc_all[grepl('PCT', colnames(RnaMetric_qc_all))] <- RnaMetric_qc_all[grepl('PCT', colnames(RnaMetric_qc_all))]/100
CombinedMetrics <- merge(RnaMetric_qc_all, MultipleMetric_qc_all, by=c('Sample','cohort',
                                                                       'PF_ALIGNED_BASES',
                                                                       'LIBRARY',
                                                                       'READ_GROUP','PF_ALIGNED_BASES'),fill=T)
CombinedMetrics <- merge(CombinedMetrics, STAR_qc_all, by=c('Sample','cohort'), all = TRUE)
CombinedMetrics$SampleFull <- CombinedMetrics$Sample
CombinedMetrics[CombinedMetrics$cohort=="Braineac",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="Braineac",]$Sample, ".*(A653.*)")[, 2]
CombinedMetrics[CombinedMetrics$cohort=="GTEx",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="GTEx",]$Sample, "(.*)_.*")[, 2]
CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample <- str_match(CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample, ".*(HRA.*)")[, 2]
CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample <- gsub('-2','',CombinedMetrics[CombinedMetrics$cohort=="TargetALS",]$Sample)

CombinedMetrics <- merge(CombinedMetrics, FastQC_qc_all_merged, by=c('Sample','cohort'), all = TRUE)
######


##### read in AMP-D ######
AMP_AD_multiMetrics_files <- list.files(path=opt$input,
                                        pattern = "CombinedMetrics.csv$", 
                                        recursive = T,
                                        full.names = T)
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
all_cohorts <- rbind.fill(CombinedMetrics, AMP_AD_multimetrics_qc_all)
all_cohorts$sampleID <- 1:nrow(all_cohorts)
#######

###### set filter based on 3 measurements ######
all_cohorts$FILTER <- "NO"
all_cohorts[!is.na(all_cohorts$PCT_CODING_BASES) & all_cohorts$PCT_CODING_BASES <= 0.1,]$FILTER <- "YES"
all_cohorts[!is.na(all_cohorts$PCT_PF_READS_ALIGNED) & all_cohorts$PCT_PF_READS_ALIGNED <= 0.6,]$FILTER <- "YES"
all_cohorts[!is.na(all_cohorts$uniquely_mapped_percent) & all_cohorts$uniquely_mapped_percent <= 60,]$FILTER <- "YES"

all_cohorts[is.na(all_cohorts$SampleFull),]$SampleFull <- all_cohorts[is.na(all_cohorts$SampleFull),]$Sample
merged_metrics_files <- paste0(opt$output,"/all_cohort_STAR_RNAseqMetrics_MultipleMetrics_FastQC.txt")
write.table(all_cohorts, merged_metrics_files, quote=F, sep="\t", row.names=F)
print(paste0("Written merged metrics to ",merged_metrics_files))

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
print(table(all_cohorts$FILTER))
print(table(all_cohorts[all_cohorts$cohort=="Brainseq",]$FILTER))

##### plot function so that it can be more easily repeated
QC_plotter <- function(all_cohorts, column, plot_pca_outliers, single_cohort =F){
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
      p <- p + geom_point(alpha=0.1)+ 
             geom_point(data=all_cohorts[all_cohorts$pcaFiltered=='YES',],  aes_string('sampleID', column, colour='cohort', shape='FILTER'))
    }else{
      p <- p + geom_point(data=all_cohorts, aes_string('sampleID', column, colour='cohort', shape='FILTER'))
    }
  }
  p <- p + theme_bw(base_size=18)+
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())+  
    xlab('Samples')
  
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
    p <- p + geom_hline(yintercept=0.1, colour="red",  linetype="dashed")
  }
  if(column == "PCT_PF_READS_ALIGNED"){
    p <- p + geom_hline(yintercept=0.6, colour="red",  linetype="dashed")
  }
  if(column == "uniquely_mapped_percent"){
    p <- p + geom_hline(yintercept=60, colour="red",  linetype="dashed")
  }
  return(p)
}
#####

##### plot STAR + multimetrics together #####

# make sure that the columns used for filtering are plotted first so that they are easy to find
columns_to_plot <- colnames(select_if(all_cohorts, is.numeric))
filter_columns <-  c("PCT_CODING_BASES","PCT_PF_READS_ALIGNED","uniquely_mapped_percent")
columns_to_plot <- columns_to_plot[!columns_to_plot %in% filter_columns]
columns_to_plot <- c(filter_columns, columns_to_plot)

pdf(paste0(opt$output,'/figures/all_MultiMetrics_QC_plots.pdf'), width=8, height=8)
for(column in columns_to_plot){
  if(column == "sampleID"){
    next
  }
  print(column)
  p <- QC_plotter(all_cohorts, column, FALSE)
  ggsave(paste0(opt$output,'/figures/QC_figures_separate/',column,'.png'), width=12, height = 8, plot=p)
  
  if(column == 'PF_READS_ALIGNED' || column == 'PCT_PF_READS_ALIGNED'){
    t <- all_cohorts[all_cohorts$cohort=="Brainseq",]
    t <- t[order(t$PCT_PF_READS_ALIGNED),]
    t$sampleID <- 1:nrow(t)
    p <- QC_plotter(t, column, FALS,TRUE)
    ggsave(paste0(opt$output,'/figures/QC_figures_separate/',column,'.Brainseq.png'), width=6, height=4)
    
  }
  print(p)
}
dev.off()
#####


##### plot STAR + multimetrics again, but highlight PCA filtered samples #####
if(!is.null(opt$pcaOutliers)){
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
#####
