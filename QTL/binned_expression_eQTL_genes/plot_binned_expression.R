library("optparse")
library(ggplot2)
library(ggpubr)
# Get command line arguments 
option_list = list(
  make_option(c("-c", "--cisQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of cis-eQTLs", metavar="character"),
  make_option(c("-t", "--transQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of trans-eQTLs", metavar="character"),
  make_option(c("-e", "--expression"), type="character",
              help="File with mean expression per gene", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt$cisQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/eQTLProbesFDR0.05-ProbeLevel.CIS.txt.gz"
opt$transQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/eQTLProbesFDR0.05-ProbeLevel.TRANS.txt.gz"
opt$expression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/MERGED_RNA_MATRIX.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.MEAN_AND_SD.txt"

##### read in data ####
cis_qtl <- read.table(gzfile(opt$cisQTLs), header=T, sep='\t')
trans_qtl <- read.table(gzfile(opt$transQTLs), header=T, sep='\t')

# make sure all are FDR < 0.05
cis_qtl <- cis_qtl[cis_qtl$FDR < 0.05,]
trans_qtl <- trans_qtl[trans_qtl$FDR < 0.05,]

expression <- read.table(opt$expression, header=T, sep='\t', row.names=1)
#####

##### per expression bin calculate the proportion of cis/trans eQTLs #####
expression$expression_bin <- cut_number(expression$mean_expression, 10)

gene_per_bin <- data.frame()
for(bin in sort(unique(expression$expression_bin))){
  genes_current_bin <- rownames(expression[expression$expression_bin==bin,])
  n_genes <- length(genes_current_bin)
  n_cis_QTL <- sum(genes_current_bin %in% cis_qtl$ProbeName)
  n_trans_QTL <- sum(genes_current_bin %in% trans_qtl$ProbeName)
  
  df <- data.frame('n_genes'=c((n_cis_QTL/n_genes)*100, 
                               (n_trans_QTL/n_genes)*100, 
                               ((n_genes-n_cis_QTL)/n_genes)*100, 
                               ((n_genes-n_trans_QTL)/n_genes)*100),
                   'bin'=c(bin, bin, bin, bin),
                   'qtl_type'=c('cis','trans','cis','trans'),
                   'genes_are_QTL'=c('yes','yes','no','no'))
  gene_per_bin <- rbind(gene_per_bin, df)
}

#####
 

#   n_genes              bin qtl_type genes_are_QTL
#1      126 [0.00777,0.0355]      cis           yes
#2       55 [0.00777,0.0355]    trans           yes

#### make the histogram ####
ggplot(gene_per_bin, aes(bin, n_genes, fill=genes_are_QTL))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=c('grey95','lightblue'))+
  theme_pubr(base_size=18)+ 
  guides(fill=FALSE)+
  facet_wrap(~qtl_type)+ 
  xlab('Expression levels')+
  ylab('Proprotion QTLs')+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
  scale_x_discrete(breaks = c('[0.00777,0.0355]','(9.43,19.3]'),
                     labels=c('low','high'))+
  theme(axis.text= element_text(colour='grey70'))
ggsave('figures/proortion_of_QTL_per_bin.pdf',width=8, height=5)  
####

#### make cummulative proportion plot ####

####