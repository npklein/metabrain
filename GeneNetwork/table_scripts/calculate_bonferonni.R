library(optparse)


### The .sh script that calls this R script has this:
# for f in ../8_GeneNetwork_predictions/*AUC.txt; do
#    echo $f;
#    Rscript calculate_bonferonni.R -i $f -o $(basename ${f%.txt}.bonferonni.txt) -t $(basename ${f%.txt}).bonSigTerms.txt;
# done

option_list <- list( 
    make_option(c("-i", "--infile"), 
        help="Infile with columns term, geneCount, pValue, auc (in this order)"),
    make_option(c("-t", "--terms"), 
        help="Termsfile with website and description of the terms"),
    make_option(c("-o", "--outfile"), 
        help="Outfile with bonf. p-value column")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(paste("Read",opt$infile))
data = read.table(opt$infile, sep ='\t', header=F, quote="")

# Here it is called term, but for GO it is an ID. However, want to keep column names the same also for HPO and kegg, which are terms
colnames(data) <- c('term', 'geneCount', 'pValue', 'auc')

data$bonferonni <- data$pValue * nrow(data)

write.table(data, file = opt$outfile, sep ='\t', quote=F, row.names=F, col.names=F)


data_bonf_sign <- data[data$bonferonni < 0.05,]
print(paste("Read",opt$terms))
terms <- read.table(opt$terms, sep='\t', comment.char="", quote="")


terms <- terms[terms$V1 %in% data_bonf_sign$term,]
write.table(terms, file=gsub('bonferonni.txt', 'bonSigTerms.txt',opt$outfile), sep='\t', col.names=F, quote=F, row.names=F)

print(paste("output written to",gsub('bonferonni.txt', 'bonSigTerms.txt',opt$outfile)))
