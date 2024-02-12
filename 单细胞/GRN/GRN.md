## GRN 基因调控网络

```properties
## 背景
http://www.seqyuan.com/GRN.html
## 数据securt多样本数据
外源刺激 -> TF(GRN) -> 转录组的变化 -> 通路/生物学过程的变化
```

#### 1.pySCENIC

```R
## 准备pySCENIC的输入文件:
## 1. 基因表达矩阵, meta cell matrix (For grnboost, step1)
## 2. TF list, download from
###   (1) 自定义
###   (2) SCENIC: https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
###   (3) AnimalTFDB: http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/TF_list_final/Homo_sapiens_TF
## 3. 参考文件:
###   (1) motif ranking tracks (.feather, 二进制文件, rows are motifs, columns are genes, values are ranks)
###   (2) motif annotation     (.tbl, txt文件, motif -> TFs)


```