

library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(viridis)
library(optparse)

option_list = list(
  make_option(c("-m", "--matrix"), type='character',
              help="Matrix to make pairwise hexbins of"),
  make_option(c("-o", "--outfile"), type='character',
              help="Outfile name"),
)
opt = parse_args(OptionParser(option_list=option_list))




hexbin_matrix <- fread(opt$matrix, header=T, na='-')

grey_grob <- rectGrob(gp=gpar(fill="white", colour='white'))
grobs <- list(grey_grob)
grobs_topSNP <- list(grey_grob)
i <- 2
columns <-colnames(hexbin_matrix)
for(colname in columns){
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
    p <- ggplot(hexbin_matrix, aes_string(colname1, colname2))+
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
   
    i <- i + 1
  }
}

opt <- list()
opt$outfile <- 'test.png'
ggsave(file = opt$outfile, arrangeGrob(grobs = grobs, ncol = length(columns)+1), width=35, height=35)  ## save plot


