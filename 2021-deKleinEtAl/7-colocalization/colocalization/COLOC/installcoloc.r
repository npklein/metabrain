if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snpStats")
install.packages("devtools", repos="https://lib.ugent.be/CRAN/")
library(devtools)
install_github("chr1swallace/coloc")