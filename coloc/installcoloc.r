if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snpStats")
install.packages("devtools")
library(devtools)
install_github("chr1swallace/coloc")