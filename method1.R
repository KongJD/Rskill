rm(list = ls())

# 数据下载 http://xena.ucsc.edu/
library(data.table)
library(dplyr)
library(tibble)

dd <- fread("./data/tcga_RSEM_gene_tpm", data.table = F)
clion <- fread("./resource/Survival_SupplementalTable_S1_20171025_xena_sp", data.table = F)

colnames(clion)[2:3] <- c("TCGA_id", "type")
index <- clion$sample %in% colnames(dd)
clion <- clion[index, 1:3]

dd <- dd[, c("sample", clion$sample)]
colnames(dd)[1] <- "gene_id"

gtf1 <- rtracklayer::import('./resource/gencode.v23.annotation.gtf')
gtf1 <- as.data.frame(gtf1)

mrna_expert <- gtf1 %>%
  filter(type == "gene", gene_type == "protein_coding") %>%
  dplyr::select(c(gene_name, gene_id)) %>%
  inner_join(dd, by = "gene_id") %>%
  dplyr::select(-gene_id) %>%
  distinct(gene_name, .keep_all = T) %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  as.data.frame() %>%
  #行名变成列名
  rownames_to_column("sample")

# table(substring(mrna_expert$sample,14,15))

mRNA_exprSet <- cbind(clion, subtype = substring(mrna_expert$sample, 14, 15), mrna_expert[, -1])

tcga_mRNA_exprSet <- mRNA_exprSet %>%
  filter(subtype != "11") %>%
  distinct(sample, .keep_all = T) %>%
  as.data.frame()

# save(tcga_mRNA_exprSet, file = "output/tcga_mRNA_exprSet.Rdata")

splitdata <- split(tcga_mRNA_exprSet, tcga_mRNA_exprSet$type)

getpancordata <- function(gene1, gene2, splitdata, method = "spearman") {
  do.call(rbind, lapply(splitdata, function(x) {
    dd <- cor.test(as.numeric(x[, gene1]), as.numeric(x[, gene2]), method = method)
    data.frame(type = x$type[1], cor = dd$estimate, p.value = dd$p.value, method = method)
  }))
}

plotdf <- getpancordata("WTAP", "SETD2", splitdata)