# RNA-seq    

### 背景知识 ###
基因表达量计算的本质：目标基因表达量相对参照系表达量的数值。        
参照的本质:        
1. 假设样本间参照的信号值应该是相同的,
2. 将样本间参照的观测值校正到同一水平;
3. 从参照的数值，校正并推算出其他观测量的值。

例如:        
- qPCR：目标基因表达量(循环数)相对看家基因表达量(循环数);        
- RNA-seq：目标基因的表达量(测序reads数)，相对样本RNA总表达量(总测序量的reads数)，这是最常用的标准。        

归一化的原因及处理原则:    
1. 基因长度
2. 测序量
3. 样本特异性(例如，细胞mRNA总量，污染等)前两者使用普通的RPKM算法就可以良好解决，关键是第三个问题，涉及到不同的算法处理。

基因表达量归一化：在高通量测序过程中，样品间在数据总量、基因长度、基因数目、高表达基因分布甚至同一个基因的不同转录本分布上存在差别。因此不能直接比较表达量，必须将数据进行归一化处理。    

RNA-seq差异表达分析的一般原则：
1. 不同样品的基因总表达量相似
2. 上调差异表达与下调差异表达整体数量相似(上下调差异平衡)
3. 在两组样品中不受处理效应影响的基因，表达量应该是相近的(差异不显著)
4. 看家基因可作为表达量评价依据(待定)

---
- DESeq2的median-of-ratios标准化方法：它假设大多数基因的表达是稳定的，并用这些“稳定基因”的表达水平来计算每个样本的标准化因子（Size Factor），以校正技术偏差（如测序深度差异）    

- 关于标准化因子的理解，可以更精确地表述：理想情况下，如果所有样本之间没有技术偏差，标准化因子会在1附近窄幅波动。如果某个样本的标准化因子显著偏离1（比如A1组因子约2.0，A2组因子约0.6），这主要反映了样本间存在的系统性技术偏差（例如总测序深度的巨大差异），而DESeq2的标准化过程正是为了校正这种偏差。标准化因子的主要作用是“纠偏”，为后续寻找真实的生物学差异奠定基础。    

- DESeq2在成功校正技术变异后，会利用广义线性模型，并基于负二项分布对计数数据进行拟合和假设检验（如Wald检验或似然比检验），来判断每个基因在组间表达量的差异是否具有统计学意义。最终，那些在经过多重检验校正后仍显示出显著差异的基因，才被认为是真正的差异表达基因。    

---
参考：
[RNA-seq中的基因表达量计算和表达差异分析](https://blog.csdn.net/weixin_30670151/article/details/96585937)        

---
步骤1-6为服务器端的上游分析，获取rawcounts，7-8为本地的下游分析，包括差异表达分析，富集分析等。
## 1.激活环境     
环境以及基因组索引是由以下链接构建的：[构建RNA-seq分析所需环境以及索引](https://github.com/y741269430/RNA-seq-for-linux/blob/main/%E6%9E%84%E5%BB%BARNA-seq%E5%88%86%E6%9E%90%E6%89%80%E9%9C%80%E7%8E%AF%E5%A2%83%E4%BB%A5%E5%8F%8A%E7%B4%A2%E5%BC%95.md)
```bash
conda activate rnaseq
```
创建文件夹，bam用于存放sam和bam文件，mapinfo用于存放比对率结果，rawcounts用于存放最终结果     
```bash
mkdir bam rawcounts mapinfo
```

## 2.写入样本名到filenames里面，用于批量运行（代码仅供参考）  
```bash
ls *1.clean* |cut -d "_" -f 1 > filenames
```

## 3.比对到mm39    
写入以下脚本到rna1_hisat2.sh中

```bash
vim rna1_hisat2.sh
```
```bash
#!/bin/bash
## mapping (hisat2) ##

mm39="/home/jjyang/downloads/genome/mm_v112/hisat2_idx/GRCm39"  # 基因组所在位置

cat filenames | while read i; 
do
nohup hisat2 -p 4 \
-x ${mm39} \
-1 ./fastq/${i}.R1.fq.gz \
-2 ./fastq/${i}.R2.fq.gz \
-S ./bam/${i}.sam 2> ./mapinfo/${i}_map.txt &

# 以下是单端比对的代码
# single end
# nohup hisat2 -p 4 \
# -x ${mm39} \
# -U ${i}_1.fq.gz \
# -S ./bam/${i}.sam 2> ./mapinfo/${i}_map.txt &
    
done
```
运行
```bash
bash rna1_hisat2.sh
```

## 4.将sam文件转换成bam文件   
写入以下脚本到rna2_sam2bam.sh中

```bash
vim rna2_sam2bam.sh
```
```bash
#!/bin/bash
## sam to bam (samtools) ##
## sorted by read name (samtools) for rowcounts ##
## sorted by position (samtools) for rmats ##

cat filenames | while read i; 
do
nohup samtools view -@ 4 -S ./bam/${i}.sam -b | samtools sort -@ 4 -n -o ./bam/${i}-sorted-name.bam &  

done
```
运行
```bash
bash rna2_sam2bam.sh
```

## 5.利用htseq-count对bam文件进行定量计算count矩阵     
写入以下脚本到rna3_htcounts.sh中   
```bash
vim rna3_htcounts.sh
```
```bash
#!/bin/bash
## calculate rawcounts (htseq-count) ##

gtf="/home/jjyang/downloads/genome/mm_v112/Mus_musculus.GRCm39.112.chr.gtf"   # gtf所在位置

cat filenames | while read i; 
do
nohup htseq-count -n 10 \
-f bam \
-r name \
-s no ./bam/${i}-sorted-name.bam ${gtf} > ./rawcounts/${i}.count &
done
```
运行
```bash
bash rna3_htcounts.sh
```

## 6.删除一些count矩阵中冗余的行 rows     
写入以下脚本到rna4_rmcounts.sh中   
```bash
vim rna4_rmcounts.sh
```
```bash
#!/bin/bash
## remove redundant rows ##

tail -n 5 ./rawcounts/* > ./rawcounts/total.info

cat filenames | while read i; 
do
sed -i '/process/d;/__/d;/retrieve/d' ./rawcounts/${i}.count & 
done
```
运行
```bash
bash rna4_rmcounts.sh
```

## hisat2/bowtie2 整合输出比对率         
输出比对率文件
```bash
conda activate rnaseq
python MappingRateOutput.py mapinfo/ alignment_res
```

- [RNAseq---Hisat2 标准输出中比对率信息解读](https://blog.csdn.net/cfc424/article/details/121666525)      

>    21800552 reads; of these:  #一共读取的reads数目，共分三部分     
>      21800552 (100.00%) were paired; of these:  #第一部分是paired-end模式下比对结果一致的结果    
>        692176 (3.18%) aligned concordantly 0 times  #在比对一致的结果中，不太合理的比对。即read1和read2都能比对上但是不合理（方向不对或者片段长度不对）    
>        20143802 (92.40%) aligned concordantly exactly 1 time   #在这些比对一致结果中，合理的比对中有92.4是reads只匹配到一处    
>        964574 (4.42%) aligned concordantly >1 times  #在这些比对的结果中，合理的比对中，有4.42reads匹配到了多处       
>        ----      
>        692176 pairs aligned concordantly 0 times; of these:  #paired-end 模式下比对不合理的比对         
>          245626 (35.49%) aligned discordantly 1 time       
>        ----        
>        446550 pairs aligned 0 times concordantly or discordantly; of these: # 剩余的reads         
>          893100 mates make up the pairs; of these:             
>            620889 (69.52%) aligned 0 times          
>            223389 (25.01%) aligned exactly 1 time          
>            48822 (5.47%) aligned >1 times          
>    98.58% overall alignment rate  #有98.58的匹配率          

## rawcounts输出解读（参考）   

https://htseq.readthedocs.io/en/master/count.html#usage  

    __no_feature: reads which could not be assigned to any feature.
    __ambiguous: reads which could have been assigned to more than one feature and hence were not counted for any of these.
    __too_low_aQual: reads which were skipped due to the -a option.
    __not_aligned: reads in the SAM file without alignment.
    __alignment_not_unique: reads with more than one reported alignment. These reads are recognized from the NH optional SAM field tag.

    __no_feature              #不能对应到任何单位类型的reads数
    __ambiguous               #不能判断落在那个单位类型的reads数
    __too_low_aQual           #低于-a设定的reads mapping质量的reads数
    __not_aligned             #存在于SAM文件，但没有比对上的reads数
    __alignment_not_unique    #比对到多个位置的reads数

---
以下内容可以本地分析      
 
## 7.合并counts文件     
```r
source('RNAseq_Flow.R')

path <- "/Users/mac/Downloads/rawcounts/"
fileNames <- dir(path, pattern = ".count$") 
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header = F, sep = "\t")})

gene_id <- data[[1]][,1]

rawdata <- Reduce(cbind, lapply(data, function(x){x <- x[,2]}))
colnames(rawdata) <- names(data)

rawdata <- cbind(gene_id, rawdata)

colnames(rawdata)[-1] <- c(paste0('Ctrl_', 1:3), 
                           paste0('Treatment_', 1:3))

write.csv(rawdata, '/Users/mac/Downloads/gene_count_matrix.csv', row.names = F)
write.table(rawdata, '/Users/mac/Downloads/gene_count_matrix.txt', row.names = F, quote = F, sep = '\t')

#### 添加annotation ####
anno <- read.delim('rmdup_20251102.txt')
anno_data <- merge(anno, rawdata, 'gene_id')
write.csv(anno_data, '/Users/mac/Downloads/anno_data.csv', row.names = F)
```
## 8.差异表达分析     
```r
#### 读取 rawdata 矩阵 ####
rawdata <- read.delim('gene_count_matrix.txt', row.names = 1)

matrix <- rawdata[,-1]
rownames(matrix) <- rawdata$gene_id

#### 绘制PCA图 ####
nor <- myNormal(matrix, myGrouplist(matrix))
plot_pca <- myPCA(nor, myGrouplist(nor))
print(plot_pca)

#### 使用DESeq2进行差异表达分析 ####
deseq2_res <- myDESeq2(matrix, 3, 3)
write.xlsx(deseq2_res, 'deseq2_res.xlsx')

#### 提取差异表达基因list ####
deg <- subset(deseq2_res, c(pvalue < 0.05 & abs(log2FC)>=0.585) )
deg <- deg[order(deg$log2FC, decreasing = T), ]
write.xlsx(deg, 'deg_data.xlsx')

#### 绘制火山图 ####
plot_vol <- myVol(deseq2_res)
print(plot_vol)

#### 绘制差异表达基因前50个上调基因与下调基因热图 ####
heat_data <- rbind(head(deg, 50), tail(deg, 50))
rownames(heat_data) <- heat_data$gene_name
heat_data <- heat_data[,-c(1:9)]
plot_heat <- myHeat(heat_data, myGrouplist(heat_data), show_rownames = T)
print(plot_heat)

#### 保存图片 ####
ggsave(plot = plot_pca, 'plot_pca.pdf', height = 5, width = 5, dpi = 300)
ggsave(plot = plot_vol, 'plot_vol.pdf', height = 8, width = 8, dpi = 300)
ggsave(plot = plot_heat, 'plot_heat.pdf', height = 8, width = 5, dpi = 300)
```
## 9.富集分析    
```r
#### 富集分析 ####
source('RNAseq_Flow.R')

# 读取全部基因list
deseq2_res <- read.xlsx('deseq2_res.xlsx')

# 构建背景基因集
deseq2_res_entrezid <- bitr(deseq2_res$gene_id, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = TRUE)

# 读取差异基因list
deg <- read.xlsx('deg_data.xlsx')

deg_entrezid <- bitr(deg$gene_id, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = TRUE)

colnames(deg_entrezid)[1] <- 'gene_id'

deg_entrezid <- left_join(deg_entrezid, deg, 'gene_id')

deg2 <- list('Up_deg' = deg_entrezid[deg_entrezid$log2FC > 0, ], 
             'Down_deg' = deg_entrezid[deg_entrezid$log2FC < 0, ])

rna <- lapply(deg2, function(x){ x <- x$ENTREZID })

# GO KEGG 富集分析
BP <- clusterProfiler::compareCluster(rna, fun = "enrichGO", ont = "BP", 
                                      OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID', 
                                      universe = deseq2_res_entrezid$ENTREZID, readable = T)
CC <- clusterProfiler::compareCluster(rna, fun = "enrichGO", ont = "CC",
                                      OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID', 
                                      universe = deseq2_res_entrezid$ENTREZID, readable = T)
MF <- clusterProfiler::compareCluster(rna, fun = "enrichGO", ont = "MF",
                                      OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID', 
                                      universe = deseq2_res_entrezid$ENTREZID, readable = T)

kegg <- setReadable(clusterProfiler::compareCluster(rna, fun = "my_enrichKEGG", organism = "mmu", 
                                                    universe = deseq2_res_entrezid$ENTREZID), 'org.Mm.eg.db', 'ENTREZID')

save(deg2, rna, BP, CC, MF, kegg, file = 'deg_GO_KEGG.RData')
```
