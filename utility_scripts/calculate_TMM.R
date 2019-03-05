#!/usr/bin/Rscript

# Normalize counts to TMM

library(edgeR)
library("optparse")

# Get command line arguments 
option_list = list(
  make_option(c("-c", "--countsFile"), type="character", default=getwd(), 
              help="path to input count file (e.g. merged STAR counts)", metavar="character"),
  make_option(c("-t", "--tmmOut"), type="character", default=getwd(),
                help="outfile with TMM normalized counts matrix", metavar="character"),
  make_option(c("-o", "--cpmOut"), type="character", default=NULL,
              help="Outfile with CPM normalized count matrix", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

countsFile <- opt$countsFile
normOut <- opt$tmmOut
cpmOut <- opt$cpmOut
samplePerc <- 0.01 # args[3] # (0.01)


counts <- read.delim(countsFile,row.names=1, check.names=FALSE)

counts_features <- counts[!rownames(counts) %in% c('N_unmapped','N_noFeature','N_multimapping',
                                                                       'N_ambiguous'),]



D <- DGEList(counts=counts_features)
cpm.D <- cpm(D)
samplePerc <- as.numeric(samplePerc)

#keep <- rowSums(cpm.D > 0.5) >= ncol(count_features)*samplePerc
#write.table(keep, file="keep_genes.txt", sep=",", row.names=TRUE)
#D <- D[keep, ,keep.lib.sizes=FALSE]
d <- calcNormFactors(D)
scalar <- d$samples$lib.size*d$samples$norm.factors/exp(mean(log(d$samples$lib.size*d$samples$norm.factors)))
scal.mat <- outer(rep(1,nrow(d$counts)), scalar)
scaled.counts <- d$counts/scal.mat

write.table(scaled.counts, file = normOut, sep = "\t", row.names=TRUE, quote=FALSE, col.names = NA)
write.table(cpm.D, file = cpmOut, sep = "\t", row.names=TRUE, quote=FALSE, col.names = NA)
