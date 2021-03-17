
library(data.table)
library(optparse)


option_list <- list(
    make_option(c("-i", "--infile"),
        help="Infile with rows genes and columns PC scores"),
    make_option(c("-o", "--outfile"),
        help="Outfile with first column scaled and then row scaled PC values")
)

opt <- parse_args(OptionParser(option_list=option_list))

pca <- data.frame(fread(paste('gunzip -c',opt$infile)), row.names=1)

pca_col_scaled <- scale(pca)

pca_row_scaled <- data.frame(t(apply(pca_col_scaled, 1, scale)))
colnames(pca_row_scaled) <- colnames(pca)

write.table(pca_row_scaled, 
            opt$outfile,
            quote=F, sep='\t', row.names = T ,col.names=T)
