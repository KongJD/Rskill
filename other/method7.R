rm(list = ls())
load(file = "data/TCGA_ggplot.Rdata")

metdata <- data.frame(TCGA_id = rownames(exprSet))

metdata <- cbind(metdata, ifelse(as.numeric(substring(metdata$TCGA_id, 14, 15)) %in% seq(1, 9), "Tumor",
                                 ifelse(as.numeric(substring(metdata$TCGA_id, 14, 15)) %in% seq(10, 29), "Normal", "")))


colnames(metdata)[2] <- "type"


load(file = "data/TCGA_BRCA_exprSet_plot.Rdata")

exprSet <- exprSet[, -c(1, 2, 3)]


caluutecor <- function(target_gene = 'MT-CO1', meth = "spearman") {
  corr_data <- data.frame()
  mt <- exprSet[, target_gene]
  for (i in 1:length(colnames(exprSet))) {
    if (colnames(exprSet)[i] == target_gene) next
    gene <- exprSet[, colnames(exprSet)[i]]
    cor <- cor.test(as.numeric(mt), as.numeric(gene), method = meth)
    corr_data[i, 1] <- target_gene
    corr_data[i, 2] <- colnames(exprSet)[i]
    corr_data[i, 3] <- cor$estimate
    corr_data[i, 4] <- cor$p.value
  }
  corr_data <- corr_data[!apply(corr_data, 1, function(x) { all(is.na(x)) }),]
  return(corr_data)
}

df <- caluutecor(target_gene = "MT-CO1", meth = "spearman")
