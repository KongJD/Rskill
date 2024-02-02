## 转录组
```properties
背景资料 
# https://mp.weixin.qq.com/s/a3y46NNNO-wardO3XWwh0w
```

#### 基本流程
```properties
流程：
1. 参考基因组、基因注释文件 下载 
2. 质控，需要fastqc及multiqc等
3. 比对：hisat2
4. 计数和归一化
5. 差异分析


流程参考：https://f1000research.com/articles/4-1070/v1 
        https://f1000research.com/articles/5-1438/v1 
        
```

```shell
step1:
cat SRR_Acc_List.txt |while read id ;do (prefetch  ${id} &);done
ls /public/project/RNA/airway/sra/*  |while read id;do ( nohup fastq-dump --gzip --split-3 -O ./ ${id} & );done
```

```shell
step2:
ls *gz |xargs fastqc -t 10
multiqc ./ 
```

```shell
step3: 去除低质量的和有接头的
trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir  $fq1 $fq2 &
```

```shell
step4:

```