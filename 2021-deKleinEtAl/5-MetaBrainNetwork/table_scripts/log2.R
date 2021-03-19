library(data.table)
library(optparse)

option_list = list(
  make_option(c("-i", "--infile"),
              help="gzipped expression file"),
  make_option(c("-o", "--outfile"), 
              help="outfile (not gzipped)")
)
opt = parse_args(OptionParser(option_list=option_list))

expression <- data.frame(fread(paste0('gunzip -c ',opt$infile), check.names=F), check.names=F, row.names = 1)
expression_log <- log2(expression+1)

write.table(expression_log, file=opt$outfile, quote=F, row.names=T, col.names=NA, sep='\t')
