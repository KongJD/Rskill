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