

library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(viridis)

zscores <- fread('2020-02-26-zsore_table.txt', header=T, na='-')
zscores_topSNP <- zscores[order(zscores$qtl_FDR),]
zscores_topSNP <- zscores_topSNP[!duplicated(zscores_topSNP$gene),]
  
colnames(zscores) <- gsub('-','\\.', colnames(zscores))


zscores_subset <- zscores[,!c('gene','snp','qtl_FDR'), with=F]


columns <- colnames(zscores_subset)


grey_grob <- rectGrob(gp=gpar(fill="white", colour='white'))
grobs <- list(grey_grob)
grobs_topSNP <- list(grey_grob)
i <- 2
for(colname in columns){
  print(colname)
  grobs[[i]] <- textGrob(colname,  gp = gpar(fontface = "bold", cex = 1))
  grobs_topSNP[[i]] <- textGrob(colname,  gp = gpar(fontface = "bold", cex = 1))
  i <- i + 1
}


for(column_index in 1:(length(columns)-1)){
  print(columns[[column_index]])
  grobs[[i]] <- textGrob(columns[[column_index]],  gp = gpar(fontface = "bold", cex = 1))
  grobs_topSNP[[i]]  <- textGrob(columns[[column_index]],  gp = gpar(fontface = "bold", cex = 1))
  i <- i + 1
  
  for(column_index2 in 1:length(columns)){
    if(column_index == column_index2){
      grobs[[i]] <- textGrob(columns[[column_index]],  gp = gpar(fontface = "bold", cex = 1))
      grobs_topSNP[[i]] <- textGrob(columns[[column_index]],  gp = gpar(fontface = "bold", cex = 1))
      i <- i + 1
      next
    }
    if(column_index2 < column_index){
      grobs[[i]] <- grey_grob 
      grobs_topSNP[[i]] <- grey_grob
      i <- i + 1
      next
    }
    colname1 <- columns[column_index]
    colname2 <- columns[column_index2]
    p <- ggplot(zscores, aes_string(colname1, colname2, colour='qtl_FDR'))+
        geom_hex()+
        theme_bw(base_size=18)+
        xlab('')+
        ylab('')+
        geom_hline(yintercept=0, lty=2, colour='red')+
        geom_vline(xintercept=0, lty=2, colour='red')+
        scale_fill_viridis()+ 
      theme(legend.position = "none",
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())
    grobs[[i]] <- p
   
     p <- ggplot(zscores_topSNP, aes_string(colname1, colname2, colour='qtl_FDR'))+
      geom_point(alpha=0.2)+
      theme_bw(base_size=18)+
      xlab('')+
      ylab('')+
      geom_hline(yintercept=0, lty=2, colour='red')+
      geom_vline(xintercept=0, lty=2, colour='red')+
      scale_colour_viridis()+ 
      theme(legend.position = "none",
            axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())
    
    grobs_topSNP[[i]] <- p
    i <- i + 1
  }
}


ggsave(file = "dataset_zscore_comparison.png", arrangeGrob(grobs = grobs, ncol = length(columns)+1), width=35, height=35)  ## save plot
ggsave(file = "dataset_zscore_comparison_topSNP.png", arrangeGrob(grobs = grobs_topSNP, ncol = length(columns)+1), width=35, height=35)  ## save plot

