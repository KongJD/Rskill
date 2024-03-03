library(dplyr)
exprSet <- gtfdata %>%
  ## 和表达量的数据交叉合并，等同于merge
  dplyr::inner_join(expr_df,by ="gene_id") %>%
  ## 去掉多余列
  dplyr::select(-gene_id) %>%
  ## 以下是为了删除重复的行(有些基因名称相同)
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>%
  ## rowMena 前置
  dplyr::select(rowMean,everything()) %>%
  ## 排序
  arrange(desc(rowMean)) %>%
  ## 去重
  distinct(gene_name,.keep_all = T) %>%
  ## 删除多余列
  dplyr::select(-rowMean)


metadata <- data.frame(TCGA_id=rownames(exprSet))

for (i in 1:nrow(metadata)) {
  ## 指示
  print(i)

  ## substring的用法，这是元素获取
  num <- as.numeric(substring(metadata[i,1],14,15))

  #如果是肿瘤，就给第2列加上Tumor
  if (num %in% seq(1,9)) {
    metadata[i,2] = "Tumor"
  }

  #如果是正常组织，就给第2列加上Normal
  if (num %in% seq(10,29)) {
    metadata[i,2] = "Normal"
  }
}
## 修改名称
colnames(metadata) <- c("TCGA_id","sample")








