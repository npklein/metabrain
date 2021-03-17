# load the data.

gene <- "CLECL1"
ct <- "macrophage"

# gene <- "CYP24A1"
# ct <- "neuron"

setwd(paste0("~/Downloads/export_ieqtl_data/", gene, "/"))

data <- read.csv("data.txt.gz", sep=",", header=TRUE, row.names = 1)
info <- read.csv("info.txt.gz", sep=",", header=FALSE, row.names = 1)

data$genotype[data$genotype == -1] <- NA
data <- data[complete.cases(data), ]
for (col in colnames(data)) {
  if (!col %in% c("genotype", "expression")) {
    data[paste0("SNPx", toupper(col))] <- data$genotype * data[col]
  }
}

plot(data$genotype, data$expression, pch = 16, col = "blue",
     xlab=info["SNPName",],
     ylab=paste0(info["HGNCName",], " expression"))
abline(lm(expression~genotype, data = data), lwd=2.5) #Add a regression line

model <- lm(expression~., data = data) 
summary(model)

plot_coeffs <- function(mlr_model) {
  coeffs <- coefficients(mlr_model)
  mp <- barplot(coeffs, col="#3F97D0", xaxt='n', main="Regression Coefficients")
  lablist <- names(coeffs)
  text(mp, par("usr")[3], labels = lablist, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.6)
}
plot_coeffs(model)

library(ggplot2)

data["group"] = factor(round(data$genotype, digits=0), levels=c(0, 1, 2))

ggplot(data=data, aes(x=data[,ct], y=expression, color=group)) + 
  geom_point(aes(shape=group), size=1.5) + 
  geom_smooth(method="lm") +
  xlab(paste0(ct, " proportion")) +
  ylab(paste0(info["HGNCName",], " expression")) +
  ggtitle("Interaction eQTL")
