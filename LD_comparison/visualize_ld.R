#!/usr/bin/Rscript

library(ggplot2)
library(data.table)
library(reshape2)
library(grid)
library(gtable)

#### Set the data: file to read. This contains the replication of MetaBrain in eQTLgen at different PC cut-offs
eqtlGen_replication <- '/groups/umcg-biogen/tmp03/output/2018-12-03-FreezeOne/analysis/output/random/z-scores/eqtlgen-FDR0.05-brain-0pcs-10pcs-20pcs-30pcs-40pcs.txt'
#

# Set the fdr threshold to filter both the meta analysis data and the replication data on
input_fdr_threshold = 0.05
#####


main <- function(eqtlGen_replication, fdr_threshold){
  # Main function for program flow
  # eqtlGen_replication: The file containing replication of MetaBrain in eqtlGen
  # fdr_threshold: fdr threshold to filter both the meta analysis data and the replication data on
  
  # Read in the table with the z-scores, extract the z-scores and calculate the replication
  # values between the meta analysis the replication dataset for 10 PC removals
  # Finally, merge everything together
  merged_zscore_df <- data.frame()
  merged_replication_df <- data.frame()

  zscore_df <- fread(eqtlGen_replication)

    
    zscore_selection_and_replication <- get_zscores_do_replication(zscore_df, study, fdr_threshold)
    
    zscore <- zscore_selection_and_replication[[1]]
    zscore$tissue <- dataset[[3]]
    merged_zscore_df <- rbind(merged_zscore_df, zscore)
    
    replication <- zscore_selection_and_replication[[2]]
    merged_replication_df <- rbind(merged_replication_df, replication)
  }
  
  # Make a dataframe that has the replication information (spearman correlation, replication rate, etc)
  # per study and per PC removal. This can be added to every facet in the plot
  studies <- unique(merged_zscore_df$study)
  datasets_text <- make_dataset_text(studies, merged_replication_df)
  
  print('Start plotting')
  plot_replication(merged_zscore_df, datasets_text, studies)
}

get_zscores_do_replication <- function(study_df, study, fdr_threshold=0.05){
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
    df_zscore_subset <- study_df[c('SNPName_0','ProbeName_0','OverallZScore_0', zscore_column, 'FDR_0',fdr_column)]
    # Save the PC in a column
    df_zscore_subset$PC <- gsub('OverallZScore_','',colnames(df_zscore_subset[4]))
    # Change the Z-score column name so they are the same ofr all the PCs
    colnames(df_zscore_subset) <- c('SNPName_0','ProbeName_0','OverallZScore_x', 'OverallZScore_y','FDR_x','FDR_y','PC')
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
    z_score_dissimilarity_score <- (sum(r(opposite_effects$OverallZScore_x-opposite_effects$OverallZScore_y)^2,na.rm=T)) / (sum((df_zscore_subset_fdr$OverallZScore_x-df_zscore_subset_fdr$OverallZScore_y)^2,na.rm=T))
    z_score_similarity_score <- round(1- z_score_dissimilarity_score,3)
    replication_scores <- rbind(replication_scores, data.frame('study'=study,
                                           'PC'=unique(df_zscore_subset_fdr$PC),
                                           'corr_coeff'=corr_coeff,
                                           'replication_rate'=replication_rate,
                                           'z_score_similarity_score'=z_score_similarity_score,
                                           'n'=nrow(!is.na(df_zscore_subset_fdr)),
                                           'opposite'=nrow(!is.na(opposite_effects))))
  }
  df_zscore$PC <- (as.numeric(df_zscore$PC)-1)*10
  replication_scores$PC <- (as.numeric(replication_scores$PC)-1)*10
  df_zscore$study <- study
  return(list(df_zscore,replication_scores))
}

make_dataset_text <- function( studies, replicate_data){
  # Make a dataframe that has the replication information 
  # Studies: List of studies that is included
  # replicate_data: Datframe with the information to print (spearman correlation, z-score similarity, etc)
  
  PCs <- c(0,10,20,30,40)
  # Need to get all combinations, e.g. 10+CMC, 10+eqtlgen, 10+brainseq, 20+CMC, 20+eqtlgen etc
  combinatins_of_PCs <-  rep(PCs, length(studies))
  combination_of_studies <- rep(studies, each=length(PCs))
  # Then add these to the dataframe
  datasets_text <- data.frame(PC = combinatins_of_PCs,
                              study = combination_of_studies)
  # Now add the labels (the replication values) from the replication data.frames we made before
  # Can do this by first merging the dataframes on study and PC
  datasets_text <- merge(datasets_text, replicate_data, by=c('study','PC'))
  
  datasets_text$label <- paste0('Spearman Correlation: ', datasets_text$corr_coeff,'\n',
                                'Z-score similarity measure: ', datasets_text$z_score_similarity_score,'\n',
                                'Replication Rate: ', datasets_text$replication_rate,'\n',
                                '# shared eQTLs: ',datasets_text$n,'\n',
                                'Opposite effects:',datasets_text$opposite)
  return(datasets_text)
}

plot_replication <- function(plotting_data, replication_text, studies){
  # Plot, using facet_grid. Make sure eQTLGen is grouped separatly as it is a blood dataset
  # plotting_data: Contains the z-scores of the different cohorts at different PC removal
  # replication_text: Contains the values that have to be printed in the plot
  
  # Scatterplot
  # make sure eQTLGen is always the last column in the plot
  studies_reordered <- c(studies[studies != 'eQTLGen'], 'eQTLGen')
  plotting_data$study <- factor(plotting_data$study, levels = studies_reordered)
  p <- ggplot(plotting_data, aes(OverallZScore_x, OverallZScore_y))+
        geom_point(alpha=0.2)+
        facet_grid(PC~study,scale='free_x')+
        theme_bw(base_size=18)+
        geom_hline(yintercept=0, lty=3, colour='red')+
        geom_vline(xintercept=0, lty=3, colour='red')+
        ylab('Meta Z-score')+
        xlab('Replication Z-score')+
        geom_text(size    = 2,
                  data    = replication_text,
                  mapping = aes(x = Inf, y = -18, label = label),
                  hjust   = 1.05,
                  vjust   = 1.5)
  outfile <- 'figures/replication_rates.png'
  ggsave(outfile, width=12, height=15)
  print(paste('Saved plot to',outfile))
  
  # Same plot but density (hex bins)
  p <- ggplot(plotting_data, aes(OverallZScore_x, OverallZScore_y))+
        geom_hex()+
        facet_grid(PC~study,scale='free_x')+
        theme_bw(base_size=18)+
        geom_hline(yintercept=0, lty=3, colour='red')+
        geom_vline(xintercept=0, lty=3, colour='red')+
        ylab('Meta Z-score')+
        xlab('Replication Z-score')+
        geom_text(size    = 2,
                  data    = replication_text,
                  mapping = aes(x = Inf, y = -18, label = label),
                  hjust   = 1.05,
                  vjust   = 1.5)+ 
        scale_fill_viridis_c(trans = "log",name = "log(#points in hex)")
  
  outfile <- 'figures/replication_rates_hexbin.png'
  ggsave(outfile, width=12, height=15)
  print(paste('Saved plot to',outfile))
  
  
  # Same but only for 10PC removal
  studies_reordered <- c(studies[studies != 'eQTLGen'], 'eQTLGen')
  plotting_data$study <- factor(plotting_data$study, levels = studies_reordered)
  p <- ggplot(plotting_data[plotting_data$PC==10,], aes(OverallZScore_x, OverallZScore_y))+
    geom_point(alpha=0.2)+
    theme_bw(base_size=50)+
    facet_grid(~study,scale='free_x')+
    geom_hline(yintercept=0, lty=3, colour='red')+
    geom_vline(xintercept=0, lty=3, colour='red')+
    ylab('Meta Z-score')+
    xlab('Replication Z-score')+
    theme(strip.text.x = element_text(size = 50))
    geom_text(size    = 6,
              data    = replication_text,
              mapping = aes(x = Inf, y = -18, label = label),
              hjust   = 1.05,
              vjust   = 1.5)
  outfile <- 'figures/replication_rates_PC10.png'
  ggsave(outfile, width=20, height=12)
  print(paste('Saved plot to',outfile))
  
  # Same plot but density (hex bins)
  p <- ggplot(plotting_data[plotting_data$PC==10,], aes(OverallZScore_x, OverallZScore_y))+
    geom_hex(aes(fill=stat(log10(count))))+
    facet_grid(~study,scale='free_x')+
    theme_pubr(base_size=50)+
    geom_hline(yintercept=0, lty=3, colour='red',size=4)+
    geom_vline(xintercept=0, lty=3, colour='red',size=4)+
    ylab('Meta Z-score')+
    xlab('Replication Z-score')+
#    geom_text(size    = 6,
#              data    = replication_text,
#              mapping = aes(x = Inf, y = -18, label = label),
#              hjust   = 1.05,
#              vjust   = 1.5)+ 
    scale_fill_viridis_c()+
    labs(fill="log10(# points in hex)")+ 
    theme(strip.text.x = element_text(size = 50))#+
    #theme(panel.background = element_rect(fill = 'black', colour = 'black'))
  
  outfile <- 'figures/replication_rates_hexbin_PC10.png'
  ggsave(outfile, width=30, height=12)
  print(paste('Saved plot to',outfile))
}

# runs only when script is run by itself, similar to Python's if __name__ == __main__
if (sys.nframe() == 0){
  main(eqtlGen_replication, input_fdr_threshold)
}

