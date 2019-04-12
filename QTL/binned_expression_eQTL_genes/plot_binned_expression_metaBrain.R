library("optparse")
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(plyr)
library(data.table)
# Get command line arguments 
option_list = list(
  make_option(c("-c", "--cisQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of cis-eQTLs", metavar="character"),
  make_option(c("-t", "--transQTLs"), type="character",
              help="eQTLProbesFDR0.05-ProbeLevel.txt.gz of trans-eQTLs", metavar="character"),
  make_option(c("-e", "--expression"), type="character",
              help="File with mean expression per gene", metavar="character"),
  make_option(c("-p", "--proteinCodingGenes"), type="character",
              help="File with protein coding genes", metavar="character"),
  make_option(c("-i", "--pLI"), type="character",
              help="File with pLI scores", metavar="character"),
  make_option(c("-o", "--out_path"), type="character",
              help="Out path for plots", metavar="character"),
  make_option(c("-z", "--ensemble_IDs"), type="character",
              help="Path to file with ensembl IDS", metavar="character"),
  make_option(c("-b", "--biosExpression"), type="character",
              help="Path to file expression data from BIOS", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- list()
opt$cisQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/eQTLProbesFDR0.05-ProbeLevel.CIS.txt.gz"
opt$transQTLs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/eQTLProbesFDR0.05-ProbeLevel.TRANS.txt.gz"
opt$expression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/MERGED_RNA_MATRIX.SampleSelection.ProbesWithZeroVarianceRemoved.QuantileNormalized.Log2Transformed.MEAN_AND_SD.txt"
opt$proteinCodingGenes <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/protein_coding_genes_ensembl84.txt"
opt$pLI <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"
opt$out_path <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/figures/"
opt$ensemble_IDs <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/ensembl_transcript_gene_IDs.txt"
opt$biosExpression <- "/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/binned_expression_eQTL_genes/data/gene_read_counts_BIOS_and_LLD_passQC.tsv.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.txt.gz"

##### read in data ####
cis_qtl <- read.table(gzfile(opt$cisQTLs), header=T, sep='\t')
trans_qtl <- read.table(gzfile(opt$transQTLs), header=T, sep='\t')
protein_coding_genes <- read.table(opt$proteinCodingGenes, sep='\t', header=T)
pLI <- read.table(opt$pLI, sep='\t', header=T)
ensemble_IDs <- fread(opt$ensemble_IDs)
BIOS_expression <- data.frame(fread(opt$biosExpression))
expression <- read.table(opt$expression, header=T, sep='\t', row.names=1)


# make sure all are FDR < 0.05
cis_qtl <- cis_qtl[cis_qtl$FDR < 0.05,]
trans_qtl <- trans_qtl[trans_qtl$FDR < 0.05,]

# genes still have version numbers, remove
rownames(expression) <- sapply(strsplit(rownames(expression),"\\."), `[`, 1)
cis_qtl$ProbeName <- sapply(strsplit(as.character(cis_qtl$ProbeName),"\\."), `[`, 1)
trans_qtl$ProbeName <- sapply(strsplit(as.character(trans_qtl$ProbeName),"\\."), `[`, 1)

# select only protein coding genes
expression_protCoding <- expression[rownames(expression) %in% protein_coding_genes$Ensembl.Gene.ID,]
cis_qtl_protCoding <- cis_qtl[cis_qtl$ProbeName %in% protein_coding_genes$Ensembl.Gene.ID,]
trans_qtl_protCoding <- trans_qtl[trans_qtl$ProbeName %in% protein_coding_genes$Ensembl.Gene.ID,]
#####

##### make expression bins #####
expression$expression_bin <- cut_number(expression$mean_expression, 10)
expression_protCoding$expression_bin <- cut_number(expression_protCoding$mean_expression, 10)
expression_protCoding$expression_bin_order = mapvalues(expression_protCoding$expression_bin, 
                                                        c("[0.00778,0.774]", "(0.774,2.86]", "(2.86,4.89]","(4.89,6.69]","(6.69,7.86]",
                                                          "(7.86,8.7]","(8.7,9.39]","(9.39,10.1]","(10.1,10.9]","(10.9,19.3]"), c(1, 2, 3,4,5,6,7,8,9,10))

#### For BIOS, calcualte mean expression and SD per gene ####
rownames(BIOS_expression) <- BIOS_expression$V1
BIOS_expression$V1 <- NULL
BIOS_expression_matrix <- BIOS_expression
BIOS_expression_matrix[BIOS_expression_matrix < 0] <- 0
BIOS_expression_matrix <- data.frame(BIOS_expression_matrix)
rownames(BIOS_expression_matrix) <- rownames(BIOS_expression)



BIOS_expression_mean <- rowMeans(BIOS_expression_matrix)
BIOS_expression_sd <- rowSds(as.matrix(BIOS_expression_matrix)) # same result as before.
BIOS_expression_mean_and_sd <- data.frame(BIOS_expression_mean, BIOS_expression_sd)
colnames(BIOS_expression_mean_and_sd) <- c('mean_expression','SD')
BIOS_expression_mean_and_sd$expression_bin <- cut_number(BIOS_expression_mean_and_sd$mean_expression, 10)
BIOS_expression_mean_and_sd$expression_bin_order = mapvalues(BIOS_expression_mean_and_sd$expression_bin, 
                                                             c("[0,0.00049]","(0.00049,0.0014]","(0.0014,0.00365]","(0.00365,0.0122]","(0.0122,0.0547]","(0.0547,0.348]",
                                                               "(0.348,2.1]","(2.1,5.71]","(5.71,8.59]", "(8.59,17.2]" ),
                                                             c(1, 2, 3,4,5,6,7,8,9,10))
####
get_gene_per_bin <- function(expression){
  # per expression bin calculate the proportion of cis/trans eQTLs and get info such as pLI and length per bin
  gene_per_bin <- data.frame()
  gene_lengths_per_bin <- data.frame()
  for(bin in sort(unique(expression$expression_bin))){
    expression_current_bin <- expression[expression$expression_bin==bin,]
    genes_current_bin <- rownames(expression_current_bin)
    n_genes <- length(genes_current_bin)
    n_cis_QTL <- sum(genes_current_bin %in% cis_qtl$ProbeName)
    n_trans_QTL <- sum(genes_current_bin %in% trans_qtl$ProbeName)
  
    mean_sd_cis <- round(mean(expression_current_bin[rownames(expression_current_bin) %in% cis_qtl$ProbeName,]$SD),2)
    mean_sd_not_cis <- round(mean(expression_current_bin[!rownames(expression_current_bin) %in% cis_qtl$ProbeName,]$SD),2)
  
    gene_lengths <- protein_coding_genes[match(rownames(expression_current_bin), protein_coding_genes$Ensembl.Gene.ID),]
    is_cis_qtl <- gene_lengths$Ensembl.Gene.ID %in% cis_qtl$ProbeName
    is_trans_qtl <- gene_lengths$Ensembl.Gene.ID %in% trans_qtl$ProbeName
  
    pLI$transcript <- gsub('\\.[0-9]+$','',pLI$transcript)
    pLI$gene <- ensemble_IDs[match(pLI$transcript, ensemble_IDs$`Transcript stable ID`),]$`Gene stable ID`
    df <- data.frame('bin'=bin, 'length'=gene_lengths$Transcript.length..including.UTRs.and.CDS.,
                   'is_cis_qtl'=is_cis_qtl,'is_trans_qtl'=is_trans_qtl,'gene'=gene_lengths$Ensembl.Gene.ID,
                   'expression'=expression_current_bin[gene_lengths$Ensembl.Gene.ID,]$mean_expression,
                   'sd'=expression_current_bin[gene_lengths$Ensembl.Gene.ID,]$mean_expression,
                   'pLI'=pLI[match(gene_lengths$Ensembl.Gene.ID, pLI$gene),]$pLI)
    gene_lengths_per_bin <- rbind(gene_lengths_per_bin, df)
  
    df <- data.frame('n_genes'=c((n_cis_QTL/n_genes)*100, 
                               (n_trans_QTL/n_genes)*100, 
                               ((n_genes-n_cis_QTL)/n_genes)*100, 
                               ((n_genes-n_trans_QTL)/n_genes)*100),
                   'bin'=c(bin, bin, bin, bin),
                   'qtl_type'=c('cis','trans','cis','trans'),
                   'genes_are_QTL'=c('yes','yes','no','no'))
    gene_per_bin <- rbind(gene_per_bin, df)
  }
  returnList <- list("gene_per_bin" = gene_per_bin, "gene_lengths_per_bin" = gene_lengths_per_bin)
  return(returnList)
}

bin_results_protCoding<- get_gene_per_bin(expression_protCoding)
bin_results<- get_gene_per_bin(expression)

#####
 
make_plots <- function(gene_per_bin, gene_lengths_per_bin, out_dir){
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
    theme(axis.text= element_text(colour='grey70'))+
    scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
                   labels=c('low','high'))
  ggsave(paste0(out_dir,'proportion_of_QTL_per_bin.png'),width=8, height=5)  
  
  
  ggplot(gene_per_bin[gene_per_bin$qtl_type=='cis',], aes(bin, n_genes, fill=genes_are_QTL))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=c('grey95','lightblue'))+
    theme_pubr(base_size=18)+ 
    guides(fill=FALSE)+
    xlab('Average brain gene expression')+
    ylab('Proportion of genes showing cis-eQTL effect')+
    scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                       labels=c('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))+
    theme(axis.text= element_text(colour='grey70'))+
    scale_x_discrete(breaks = c('[0.00778,0.774]','(10.9,19.3]'),
                     labels=c('low','high'))
  ggsave(paste0(out_dir,'proportion_of_QTL_per_bin_only_cis.png'),width=8, height=5)  
  
  ####

  #### per bin plot the length of the genes for the QTL genes and non-QTL genes ####
  p1 <- ggplot(gene_lengths_per_bin, aes(log(length), fill=is_cis_qtl))+
    geom_histogram(position="identity", alpha=0.5)+
    facet_wrap(~bin, ncol=2)+
    theme_pubr(base_size=18)+
    scale_fill_manual(values=c('grey70','blue'))+
    labs(fill="is cis QTL")
  ggsave(paste0(out_dir,'length_per_expression_bin_cis.pdf'),width=8, height=12)  
  
  
  p2 <- ggplot(gene_lengths_per_bin, aes(log(length), fill=is_trans_qtl))+
    geom_histogram(position="identity", alpha=0.5)+
    facet_wrap(~bin, ncol=2)+
    theme_pubr(base_size=18)+
    scale_fill_manual(values=c('grey70','blue'))+
    labs(fill="is trans QTL")
  ggsave(paste0(out_dir,'length_per_expression_bin_trans.pdf'),width=8, height=12)  
  
  pdf(paste0(out_dir,'length_per_expression_bin.pdf'),width=12, height=12)  
  grid.arrange(p1, p2, nrow = 1)
  dev.off()
  ####
  
  
  #### writing genes from highest bin for GO analysis ####
  highest_bin <- gene_lengths_per_bin[gene_lengths_per_bin$bin == "(10.9,19.3]",]
  write.table(as.character(highest_bin[highest_bin$is_cis_qtl,]$gene), 'cis_qtl_genes_highest_bin.txt',
              quote=F, row.names=F)
  write.table(as.character(highest_bin[!highest_bin$is_cis_qtl,]$gene), 'not_cis_qtl_genes_highest_bin.txt',
              quote=F, row.names=F)
  ####
  
  #### plot sd per bin for cis vs non cis ####
  p1 <- ggplot(gene_lengths_per_bin, aes(sd, fill=is_cis_qtl))+
    geom_histogram(position='identity')+
    theme_pubr(base_size=18)+ 
    scale_fill_manual(values=c('grey95','lightblue'))+
    facet_wrap(~bin, ncol=2,scale='free_x')
  
  p2 <- ggplot(gene_lengths_per_bin, aes(sd, fill=is_trans_qtl))+
    geom_histogram(position='identity')+
    theme_pubr(base_size=18)+ 
    scale_fill_manual(values=c('grey95','lightblue'))+
    facet_wrap(~bin, ncol=2,scale='free_x')
  pdf(paste0(out_dir,'sd_per_expression_bin.pdf'),width=12, height=12)  
  grid.arrange(p1, p2, nrow = 1)
  dev.off()
  ####
  
  #### Make pLI plot ####
  gene_lengths_per_bin$is_cis_qtl <- factor(gene_lengths_per_bin$is_cis_qtl, levels=c(TRUE,FALSE))
  ggplot(gene_lengths_per_bin, aes(is_cis_qtl,pLI, fill=is_cis_qtl))+
    geom_violin()+
    stat_summary(fun.y = median, geom = "point", shape=23, size = 3)+
    theme_bw(base_size=18)+
    facet_grid(~bin)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top')+
    scale_fill_manual(values=c('lightblue','grey95'))+
    xlab('')
  ggsave(paste0(out_dir,'pLI_per_bin.png'),height=5, width=12)
  
  ####
}

make_plots(bin_results_protCoding$gene_per_bin, bin_results_protCoding$gene_lengths_per_bin, paste0(opt$out_path,'proteinCoding/'))
make_plots(bin_results$gene_per_bin, bin_results$gene_lengths_per_bin, paste0(opt$out_path,'all/'))

expression_protCoding$name <- 'MetaBrain'
BIOS_expression_mean_and_sd$name <- 'BIOS'
expression_metaBrain_bios <- rbind(expression_protCoding,BIOS_expression_mean_and_sd)
ggplot(expression_metaBrain_bios, aes(expression_bin_order, SD,fill=name))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  xlab('Expression bin')
ggsave(paste0(opt$out_path,'proteinCoding/sd_comparison.png'),width=8, height=8)
ggplot(expression_metaBrain_bios, aes(expression_bin_order, mean_expression,fill=name))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  xlab('Expression bin')
ggsave(paste0(opt$out_path, 'proteinCoding/mean_comparison.png'),width=8, height=8)

