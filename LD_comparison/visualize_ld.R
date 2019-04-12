#!/usr/bin/Rscript
library(optparse)
library(ggplot2)
library(data.table)
library(reshape2)
library(grid)
library(gtable)
library(ggpubr)
library('scales')


option_list = list(
  make_option(c("-e", "--eqtlGenReplication"), type="character",
              help="File containing the replication of MetaBrain in eQTLgen at different PC cut-offs", metavar="character"),

  make_option(c("-l", "--ld_scores"), type="character",
              help="File containing the LD scores between the top SNPs for each gene of metaBrain that is also found in eqtlGen (fdr < 0.05)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#### Set the data: file to read ####
# This contains the replication of MetaBrain in eQTLgen at different PC cut-offs
eqtlGen_replication <- opt$eqtlGenReplication
eqtlGen_replication <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/replications/data/eqtlgen-FDR0.05-brain-0pcs-10pcs-20pcs-30pcs-40pcs.txt'

# This contains the LD scores between the top SNPs for each gene of metaBrain that is also found in eqtlGen (fdr < 0.05)
ld_scores_file <- opt$ld_scores
ld_scores_file <- '/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/LD_comparison/metaBrain_eqtlGen_topSNP_LD.txt'

#####

# Set the fdr threshold to filter both the meta analysis data and the replication data on
fdr_threshold = 0.05
#####


main <- function(eqtlGen_replication, ld_scores_file, fdr_threshold){
  # Main function for program flow
  # eqtlGen_replication: The file containing replication of MetaBrain in eqtlGen
  # fdr_threshold: fdr threshold to filter both the meta analysis data and the replication data on
  
  # Read in the table with the z-scores, extract the z-scores and calculate the replication
  # values between the meta analysis the replication dataset for 10 PC removals
  # Finally, merge everything together
  merged_zscore_df <- data.frame()
  merged_replication_df <- data.frame()

  zscore_df <- fread(eqtlGen_replication)

  # get the zscore selection and replication for the difference PC corrected data, can select PC10 later
  zscore_selection_and_replication <- get_zscores_do_replication(zscore_df, fdr_threshold)
    
  zscores <- zscore_selection_and_replication[[1]]
  # select disconcordant/concordant effects
  zscores$concordant <- sign(zscores$OverallZScore_x) == sign(zscores$OverallZScore_y)
  
  replication <- zscore_selection_and_replication[[2]]

  plot_replication(zscores, replication)
  
  
  # Read the LD data
  ld_scores <- fread(ld_scores_file)
  # Some of the comparisons could not be made and have error instead of value, remove these
  ld_scores <- ld_scores[!grepl('error', ld_scores$`D'`),]
  
  # For the concordant and disconcordant eQTLs, plot the D' and R2 for the top SNPs of metabrain vs eQTLgen
  plot_ld(zscores, ld_scores)
}

get_zscores_do_replication <- function(study_df, fdr_threshold=0.05){
  # Select columns + change to long format, add study as column, then merge. Also since we are looping already, do the correlation calculations etc
  # study_df: table with z-scores for 1 cohort at different PC removal steps
  # study: Name of the study, e.g. Brainseq
  # fdr_threshold: Threshold to select significant hits in *both* meta analysis and replication cohort
  
  study_df <- data.frame(study_df)
  zscore_columns <- colnames(study_df)[grepl('OverallZScore',colnames(study_df))]
  # make df for each of the z-scores, then merge after. Also calculate the correlations etc
  df_zscore <- data.frame()
  replication_scores <- data.frame()
  for(zscore_column in zscore_columns){
    # comparing overallZscore_0 to everything so don't need to add it here
    if(zscore_column=='OverallZScore_0'){
      next
    }
    
    fdr_column <- gsub('OverallZScore','FDR', zscore_column)
    alleleAssed_column <- gsub('OverallZScore','AlleleAssessed', zscore_column)
    snpType_column <- gsub('OverallZScore','SNPType', zscore_column)
    df_zscore_subset <- study_df[c('SNPName_0','ProbeName_0','OverallZScore_0', zscore_column, 
                                   'FDR_0',fdr_column,'AlleleAssessed_0',alleleAssed_column,
                                   'SNPType_0',snpType_column)]
    # Save the PC in a column
    df_zscore_subset$PC <- gsub('OverallZScore_','',colnames(df_zscore_subset[4]))
    # Change the Z-score column name so they are the same ofr all the PCs
    colnames(df_zscore_subset) <- c('SNPName_0','ProbeName_0','OverallZScore_x', 'OverallZScore_y','FDR_x','FDR_y',
                                    'AlleleAssessed_x','AlleleAssessed_y','SNPType_x','SNPType_y','PC')
    df_zscore_subset_fdr <- df_zscore_subset[df_zscore_subset$FDR_x <= fdr_threshold & df_zscore_subset$FDR_y <= fdr_threshold,]
    df_zscore_subset_fdr <- df_zscore_subset_fdr[!is.na(df_zscore_subset_fdr$PC),]
    df_zscore <- rbind(df_zscore, df_zscore_subset_fdr)
    
    # Calculate the correlation
    corr_coeff <-round(cor(df_zscore_subset_fdr$OverallZScore_x, df_zscore_subset_fdr$OverallZScore_y, method='spearman',use='complete.obs'),3)
    
    # Calculate the replication rate
    opposite_effects <- df_zscore_subset_fdr[df_zscore_subset_fdr$OverallZScore_x * df_zscore_subset_fdr$OverallZScore_y < 0,]
    same_effects <- df_zscore_subset_fdr[df_zscore_subset_fdr$OverallZScore_x * df_zscore_subset_fdr$OverallZScore_y > 0,]
    opposite_effects_rate <- nrow(opposite_effects[!is.na(opposite_effects$SNPName_0),])/nrow(df_zscore_subset_fdr[!is.na(df_zscore_subset_fdr$SNPName_0),])
    replication_rate <- round(1-opposite_effects_rate,3)
    
    # Calculate the dissimilarity score
    z_score_dissimilarity_score <- (sum((opposite_effects$OverallZScore_x-opposite_effects$OverallZScore_y)^2,na.rm=T)) / (sum((df_zscore_subset_fdr$OverallZScore_x-df_zscore_subset_fdr$OverallZScore_y)^2,na.rm=T))
    z_score_similarity_score <- round(1- z_score_dissimilarity_score,3)
    replication_scores <- rbind(replication_scores, data.frame('PC'=unique(df_zscore_subset_fdr$PC),
                                           'corr_coeff'=corr_coeff,
                                           'replication_rate'=replication_rate,
                                           'z_score_similarity_score'=z_score_similarity_score,
                                           'n'=nrow(!is.na(df_zscore_subset_fdr)),
                                           'opposite'=nrow(!is.na(opposite_effects))))
  }
  df_zscore$PC <- (as.numeric(df_zscore$PC)-1)*10
  replication_scores$PC <- (as.numeric(replication_scores$PC)-1)*10
  return(list(df_zscore,replication_scores))
}

plot_replication <- function(zscores, replication, PCs=10){
  print('Start plotting')
  # Plot, using facet_grid. Make sure eQTLGen is grouped separatly as it is a blood dataset
  # zscores: Contains the z-scores of the different cohorts at different PC removal
  # replication: Contains the values that have to be printed in the plot
  
  replication <- replication[replication$PC==PCs,]
  replication$label <- paste0('Replication Rate: ', replication$replication_rate,'\n',
                   '# shared eQTLs: ',replication$n,'\n',
                   'Opposite effects:',replication$opposite)
  # Scatterplot eqtlgen vs metabrain
  p <- ggplot(zscores[zscores$PC==PCs,], aes(OverallZScore_x, OverallZScore_y, colour=concordant))+
    geom_point(alpha=0.2)+
    theme_bw(base_size=18)+
    geom_hline(yintercept=0, lty=3, colour='red')+
    geom_vline(xintercept=0, lty=3, colour='red')+
    ylab('Meta Z-score')+
    xlab('Replication Z-score')+
    theme(strip.text.x = element_text(size = 50))+
    annotate("text",
             size    = 3,
             x = max(zscores[zscores$PC==PCs,]$OverallZScore_x), 
             y = min(zscores[zscores$PC==PCs,]$OverallZScore_y), 
             label = replication$label,
              hjust   = 0.95,
              vjust   = 0.1)+ 
    scale_colour_brewer(palette='Set1')+
    guides(colour=F)
  outfile <- 'figures/replication_rates_PC10.png'
  ggsave(outfile, width=6, height=6)
  print(paste('Saved plot to',outfile))
  
  
  zscores_subset <- zscores[zscores$ProbeName_0=='ENSG00000065621' | zscores$ProbeName_0=='ENSG00000134061',]

  
}

plot_ld <- function(zscores, ld_scores, PCs = 10){
  zscores <- zscores[zscores$PC==PCs,]
  # remove the version number of gene because ld data does not have it
  zscores$ProbeName_0 <- gsub('\\.[0-9]+','',zscores$ProbeName_0)
  # For the genes in ld_scores, add column if they were concordant or disconcordant direction
  zscores$SNPName_0 <- NULL
  zscores$PC <- NULL
  ld_scores_merged <- merge(ld_scores, zscores, by.x=c('gene'), by.y=c('ProbeName_0'))
  colnames(ld_scores_merged) <- c('gene','topSNP_metaBrain','topSNP_eqtlGen',"D'",'R2',
                                  'OverallZScore_MetaBrain','OverallZScore_eqtlGen',
                                  'FDR_MetaBrain','FDR_eqtlGen','AlleleAssessed_MetaBrain',
                                  'AlleleAssessed_eqtlGen','SNPType_MetaBrain','SNPType_eqtlGen',
                                  'concordant')
  hgnc <- fread('ensembl_to_hgnc.txt')
  ld_scores_merged$hgnc <- hgnc[match(ld_scores_merged$gene,hgnc$`Gene stable ID`),]$`HGNC symbol`
  # for som genes there is no concorandt info because they are > 0.05 in metabrain, remove these
  ld_scores_merged <- ld_scores_merged[!is.na(ld_scores_merged$concordant),]
  
  # because error messages were in the column they were considered characters, convert to numeric
  ld_scores_merged$`D'` <- as.numeric(ld_scores_merged$`D'`)
  ld_scores_merged$R2 <- as.numeric(ld_scores_merged$R2)
  
  # make the plot
  p <- ggplot(ld_scores_merged, aes(x=`D'`, y=R2, colour=concordant))+
    geom_point(alpha=0.5)+
    xlab("D'")+
    theme_pubr(base_size=18)+
    scale_colour_brewer(palette='Set1')
  
  
  pdf('figures/ld_comparison_density.pdf',width=12, height=8)
  ggExtra::ggMarginal(p, type = "density",  groupFill = TRUE)
  dev.off()
  
  pdf('figures/ld_comparison_histogran.pdf',width=12, height=8)
  ggExtra::ggMarginal(p, type = "histogram",  groupFill = TRUE)
  dev.off()
  
  write.table(ld_scores_merged, file='LD_scores_topSNP_metaBrain_eqtlGen.txt',sep='\t','quote'=F, row.names = F)

  ld_scores_merged_disconcordant <- ld_scores_merged[ld_scores_merged$concordant==F & ld_scores_merged$R2 == 1,]
  write.table(ld_scores_merged_disconcordant, file='LD_scores_topSNP_metaBrain_eqtlGen.disconcordant.txt',sep='\t','quote'=F, row.names = F)
  
  # merge zscores and ld
  ld_zscore <- merge(zscores, ld_scores_merged, by.x='ProbeName_0', by.y='gene')
  
  # Look at the zscore
  p <- ggplot(ld_zscore, aes(x=abs(OverallZScore_x), y=R2, colour=concordant.x))+
    geom_point()+
    xlab('metabrain_zscore')+
    theme_bw(base_size=18)
  pdf('figures/ld_metbrain_zscore.pdf',width=12, height=8)
  ggExtra::ggMarginal(p, type = "density",  groupFill = TRUE)
  dev.off()
  
  
  # also look at tss distance
  eqtlGen_distance <- fread('/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/TSS_distance/tssdist-eqtlgen.txt')
  colnames(eqtlGen_distance) <- c('pval','gene','snp','eqtlgen_snp_distance')
  ld_zscore <- merge(ld_zscore, eqtlGen_distance, by.x='ProbeName_0',by.y='gene')
  
  metaBrain_distance <- fread('/Users/NPK/UMCG/projects/biogen/cohorts/joined_analysis/TSS_distance/eQTLsFDR0.05-TSS-10pcs.topfxpergene.txt')
  colnames(metaBrain_distance)[4] <- 'metaBrain_snp_distance'
  metaBrain_distance$ProbeName <- gsub('\\.[0-9]+','', metaBrain_distance$ProbeName)
  ld_zscore <- merge(ld_zscore, metaBrain_distance, by.x='ProbeName_0', by.y='ProbeName')
  
  ld_zscore$dist_diff <- abs(abs(ld_zscore$eqtlgen_snp_distance)-abs(ld_zscore$metaBrain_snp_distance))
  
  p <- ggplot(ld_zscore, aes(x=dist_diff, y=R2,colour=concordant.x))+
    geom_point(alpha=0.2)+
    theme_bw(base_size=18)+
    xlab('absolute difference in distance snp-gene between eqtlGen and metaBrain')
  pdf('figures/ld_metbrain_dist_diff.pdf',width=12, height=8)
  ggExtra::ggMarginal(p, type = "density",  groupFill = TRUE)
  dev.off()
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main(eqtlGen_replication, ld_scores_file, fdr_threshold)
}

