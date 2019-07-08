library(data.table)

input_file <- '../test_less_genes/output/5_covariatesRemoved/1000_genes_test.duplicateSamplesRemoved_extractedColumnsnoVarianceRowsRemoved.DESeqNorm.CovariatesRemoved.txt'
expression <- fread(input_file,
                    check.names=F)
expression <- data.frame(expression, check.names=F)
rownames(expression) <- expression$`-`
expression$`-` <- NULL


tissue_info <- fread('brain.phenotypes.txt')

calc_coexpr <- function(expr,expr_t,  gene, tissue=''){

  if(nchar(tissue) > 0){
    print(paste('subset on',tissue))
    tissue_samples <- tissue_info[which(tissue_info$BroadBrainRegion==tissue),]
    expr <- expr[colnames(expr) %in% tissue_samples$rnaseq_id]
    expr_t <- expr_t[rownames(expr_t) %in% tissue_samples$rnaseq_id,]
  }
  gene_expr <- as.numeric(expr[gene,])

  print('Calculate correlation')
#  coexpress <- apply(expr, 1, function(x){
#    cor.test(gene_expr, x, method="spearman")
#e  })

  result <- do.call(rbind, 
                    lapply(1:nrow(expr), function(x){
                        cor.result <- cor.test(gene_expr, as.numeric(expr[x,]))
                        pvalue <- cor.result$p.value
                        estimate <- cor.result$estimate
                        df <- data.frame(pvalue=pvalue, estimate=estimate)
                        rownames(df) <- rownames(expr)
                        return(df)}))

  coexpress <- cor.test(gene_expr, expr_t, method="spearman")

  outname <- paste0(gene,'_spearman_',ncol(expr),'_samples')
  if(nchar(tissue) > 0){outname <- paste0(outname,'_',tissue)}
  outfile <- paste0('co_expression_files/',outname,'.txt')
  write.table(t(coexpress), outfile,
              quote=F, col.names=F, sep='\t')
  print(paste('Written to', outfile))

}

expression_trans <- t(expression)
for(gene in c('ENSG00000069329')){
  print(paste('start',gene))
  calc_coexpr(expression, expression_trans,gene)
  for(tissue in c("Basal ganglia", "Cerebellum", "Cortex", "Limbic system","Spinal cord")){
    alc_coexpr(expression, expression_trans,gene, tissue = tissue)
  }
}
