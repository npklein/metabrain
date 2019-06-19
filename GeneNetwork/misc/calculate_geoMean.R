#!/usr/bin/env Rscript
library("optparse")
library(data.table)
library(EnvStats)
option_list = list(
  make_option(c("-e", "--expr"), type="character",
              help="Expression data (gzipped only)", metavar="character"),
    make_option(c("-o", "--out"), type="character", 
              help="output file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


expression <- data.frame(fread(paste0('gunzip -c ', opt$expr)))
rownames(expression) <- expression$Rowname
expression$Rowname <- NULL
expression <- expression+1

geoMean <- data.frame(apply(expression, 1, geoMean))
today <- as.character(Sys.Date())
geoMean[,today] <- rownames(geoMean)

colnames(geoMean)[1] <- 'Col1'
geoMean <- geoMean[,c(today, 'Col1')]
geoMean$Col1 <- geoMean$Col1 -1
write.table(geoMean, file = opt$out,quote=F,row.names=F, sep='\t')

