
library(data.table)
library(ggplot2)
library("dplyr")
library(stringr)
library("optparse")
library(RColorBrewer)
library(tidyr)
library("readxl")
library(tidyverse)



# Get command line arguments 
option_list = list(
  make_option(c("-i", "--input"), type="character", default=getwd(), 
             help="path to input file (phenotype table)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=getwd(),
                help="path to output dir", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=getwd(),
              help="File with samples to include", metavar="character"),
  make_option(c("-p", "--populations"), type="character", default=getwd(),
              help="File with population numbers", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

test <- function(){
  opt <- list()
  opt$input <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2020-02-03-phenotype-table/2020-03-09.brain.phenotypes.txt'
  opt$samples <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2019-01-14-QC/2020-02-03-all-expression-samples.txt'
  opt$output <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2020-02-03-phenotype-table/figures/'
  opt$population <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/2019-11-06-freeze2dot1/2020-02-03-phenotype-table/individuals_population.txt'
  return(opt)
}
dir.create(file.path(opt$output), showWarnings = FALSE)


populations_afterQC <- read.table(opt$population, sep="\t", header=T)
populations_afterQC$dataset <- ""
populations_afterQC[grepl("AMPAD", populations_afterQC$Genotype.dataset),]$dataset <- "AMP-AD"
populations_afterQC[grepl("Bipseq", populations_afterQC$Genotype.dataset),]$dataset <- "PsychEncode"
populations_afterQC[grepl("LIBD", populations_afterQC$Genotype.dataset),]$dataset <- "PsychEncode"
populations_afterQC[grepl("Braineac", populations_afterQC$Genotype.dataset),]$dataset <- "Braineac"
populations_afterQC[grepl("TargetALS", populations_afterQC$Genotype.dataset),]$dataset <- "TargetALS"
populations_afterQC[grepl("UCLA_ASD", populations_afterQC$Genotype.dataset),]$dataset <- "PsychEncode"
populations_afterQC[grepl("NABEC", populations_afterQC$Genotype.dataset),]$dataset <- "NABEC"
populations_afterQC[grepl("GVEX", populations_afterQC$Genotype.dataset),]$dataset <- "PsychEncode"
populations_afterQC[grepl("GTEx", populations_afterQC$Genotype.dataset),]$dataset <- "GTEx"
populations_afterQC[grepl("CMC", populations_afterQC$Genotype.dataset),]$dataset <- "CMC"
populations_afterQC[grepl("CMC_HBCC", populations_afterQC$Genotype.dataset),]$dataset <- "PsychEncode"
populations_afterQC[grepl("ENA", populations_afterQC$Genotype.dataset),]$dataset <- "ENA"

populations_resummed <- populations_afterQC %>% group_by(dataset) %>% summarise(EUR=sum(EUR), AFR=sum(AFR), EAS=sum(EAS), SAS=sum(SAS, AMR=sum(AMR)))


print("read phenotypes")
samples <- read.table(opt$samples)
phenotypes <- fread(opt$input,sep="\t",quote="")
phenotypes <- phenotypes[phenotypes$rnaseq_id %in% samples$V1,]
phenotypes[phenotypes$predicted.brain.region == "no_predicted_region_available",]$predicted.brain.region <- "region not predicted"

phenotypes[is.na(phenotypes$Diagnosis),]$Diagnosis <- 'Not available'

colourCount = length(unique(phenotypes$Diagnosis))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
dark2ExtendedPallete <- getPalette(colourCount)
dark2ExtendedPallete[length(dark2ExtendedPallete)] <- "grey90"

phenotypes[grepl('Dementia',phenotypes$Diagnosis),]$Diagnosis <- 'Dementia'
levels <- unique(phenotypes$Diagnosis)
levels <- levels["Not available" != levels]
phenotypes$Diagnosis <- factor(phenotypes$Diagnosis, levels=c(levels,"Not available"))
phenotypes$MetaCohort <- fct_infreq(phenotypes$MetaCohort)
print("plot brain region")
ggplot(phenotypes, aes(predicted.brain.region,fill=Diagnosis,colour='yes'))+
  geom_bar()+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F)+ 
  coord_flip()
ggsave(paste0(opt$output,"/diagnosis_legend.pdf"), width=10, height=5)

ggplot(phenotypes, aes(predicted.brain.region,fill=Diagnosis,colour='yes'))+
  geom_bar()+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F)+ 
  coord_flip()+
  guides(fill=F)
ggsave(paste0(opt$output,"/diagnosis_noLegend.pdf"), width=9, height=6)

###### color ######
wong_hex <- c(rgb(0,0,0, maxColorValue = 255), 
              rgb(230,159, 0, maxColorValue = 255), 
              rgb(86,180,233, maxColorValue = 255), 
              rgb(204,121,167, maxColorValue = 255), 
              rgb(240,228,66, maxColorValue = 255), 
              rgb(0,114,178, maxColorValue = 255), 
              rgb(213, 94,0, maxColorValue = 255), 
              rgb(0,158,115, maxColorValue = 255))
wong_hex_metacohort <- c(wong_hex, "#A9A9A9")
wong_hex_population <- c(wong_hex_metacohort[6], wong_hex_metacohort[2], wong_hex_metacohort[8], wong_hex_metacohort[3])
wong_order_brain <- c(wong_hex[4], wong_hex[5],wong_hex[2],wong_hex[3],
                      wong_hex[8],wong_hex[7],wong_hex[6])

######

##### populations #####
populations_resummed_melt <- melt(populations_resummed)


print("plot populations")
populations_resummed_melt$variable <- factor(populations_resummed_melt$variable, levels=c('SAS','EAS','AFR','EUR'))
ggplot(populations_resummed_melt, aes(reorder(dataset, -value), value,fill=variable,colour='yes'))+
  geom_bar(stat='identity',width = 0.5)+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of individuals')+
  scale_fill_manual(values=wong_hex_population)+
  scale_colour_manual(values='black')+
  guides(colour=F,fill = guide_legend(reverse = TRUE))+ 
  coord_flip()+
  theme(legend.position = c(0.7, 0.8))+
  labs(fill="Population")
ggsave(paste0(opt$output,"/2020-06-25-fig1B-populations-wrong-grouping.pdf"), width=6, height=6)

#######



phenotypes <- phenotypes[!phenotypes$predicted.brain.region=="region not predicted",]

phenotypes[phenotypes$predicted.brain.region=='amygdala',]$predicted.brain.region <- 'Amygdala'
phenotypes[phenotypes$predicted.brain.region=='basalganglia',]$predicted.brain.region <- 'Basal ganglia'
phenotypes[phenotypes$predicted.brain.region=='cerebellum',]$predicted.brain.region <- 'Cerebellum'
phenotypes[phenotypes$predicted.brain.region=='cortex',]$predicted.brain.region <- 'Cortex'
phenotypes[phenotypes$predicted.brain.region=='hippocampus',]$predicted.brain.region <- 'Hippocampus'
phenotypes[phenotypes$predicted.brain.region=='hypothalamus',]$predicted.brain.region <- 'Hypothalamus'
phenotypes[phenotypes$predicted.brain.region=='spinalcord',]$predicted.brain.region <- 'Spinal cord'
phenotypes$predicted.brain.region <- factor(phenotypes$predicted.brain.region, levels=c('Amygdala','Hippocampus','Hypothalamus',
                                                                                        'Spinal cord','Basal ganglia',
                                                                                        'Cerebellum','Cortex'))
ggplot(phenotypes, aes(reorder(MetaCohort,MetaCohort,
                               function(x)-length(x)),fill=predicted.brain.region,colour='yes'))+
  geom_bar(width = 0.5)+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=wong_order_brain)+
  scale_colour_manual(values='black')+
  guides(colour=F,fill = guide_legend(reverse = TRUE))+ 
  coord_flip()+
  theme(legend.position = c(0.7, 0.75))+
  labs(fill="Brain region")
ggsave(paste0(opt$output,"/2020-06-25-fig1A-brain-regions-wrong-grouping.pdf"), width=6, height=6)

ggplot(phenotypes, aes(MetaCohort,fill=`Brain region`,colour='yes'))+
  geom_bar()+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F, fill=F)+ 
  coord_flip()
  
  ggsave(paste0(opt$output,"/region_per_cohort_noLegend.pdf"), width=6, height=9)



colourCount = length(unique(phenotypes$MetaCohort))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
dark2ExtendedPallete <- getPalette(colourCount)
ggplot(phenotypes, aes(MetaCohort,fill=MetaCohort,colour='yes'))+
  geom_bar()+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F)+ 
  coord_flip()+
  guides(fill=F)
  
ggsave(paste0(opt$output,"/samples_per_cohort.pdf"), width=6, height=9)



phenotypes <- data.frame(phenotypes)

phenotypes[phenotypes$Age=="90+" & !is.na(phenotypes$Age),]$Age <- "90"
phenotypes$Age <- as.numeric(phenotypes$Age)
colourCount = length(unique(phenotypes$`Brain region`))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
dark2ExtendedPallete <- getPalette(colourCount)
phenotypes_no_na_age <- phenotypes[!is.na(phenotypes$Age),]
ggplot(phenotypes_no_na_age, aes(Age,fill=`MetaCohort`,colour='yes'))+
  geom_histogram()+ 
  theme_bw(base_size=20)+
  xlab('')+
  ylab('Number of RNAseq samples')+
  scale_fill_manual(values=dark2ExtendedPallete)+
  scale_colour_manual(values='black')+
  guides(colour=F)
ggsave(paste0(opt$output,"/age.pdf"), width=9, height=6)
