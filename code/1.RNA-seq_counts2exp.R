##### 根据定量文件生成基因和转录本的count矩阵
setwd(r"{E:\项目-YJJ\save\}")

# 批量创建文件夹
dir.create('E:/项目-YJJ/save/Output/merged_result/', recursive = T)
dir.create('E:/项目-YJJ/save/Output/COND1/', recursive = T)

group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

base_path <- 'E:/项目-YJJ/save/Output/COND1/'

subfolders <- c('func_result', 'gsea_result', 'mRNA_result', 'diff_result')

for (g in group) {
  for (sub in subfolders) {
    folder_path <- file.path(base_path, g, sub)
    dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  }
}

#### 输入count文件
data <- read.csv('Output/gene_count_matrix.csv')


test <- read.delim('ref_transcripts_expression.txt')
test <- test[,c(6,9)]

test <- left_join(test, data)
test <- na.omit(test)

test <- test[,-2]

colnames(test)[1] <- 'transcript_id'

write.csv(test, 'Output/transcript_count_matrix.csv', row.names = F)
write.table(test, 'Output/transcript_count_matrix.txt', row.names = F, quote = F, sep = '\t')

##### 生成fpkm.txt ####

myGrouplist <- function(x){
  num <- stringr::str_split_fixed(colnames(x), "_", n = 2)[,1]
  num <- factor(num, levels = unique(num)); num
  num <- data.frame(table(num)); num
  
  grouplist <- rep(num[1:nrow(num),1], num[1:nrow(num),2]); grouplist
}
myNormal <- function(x, grouplist){
  colData <- data.frame(row.names = colnames(x), grouplist)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                        colData = colData,
                                        design = ~ grouplist)
  keep <- rowSums(DESeq2::counts(dds)) >= 0
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)  
  nor <- DESeq2::counts(dds, normalized = T)
  nor2 <- as.data.frame(nor)
}

rownames(data) <- data$gene_id
data <- data[,-1]

nor <- myNormal(data, myGrouplist(data))

write.table(nor, 'Output/fpkm.txt', row.names = T, quote = F, sep = '\t')


#####
matrix <- read.delim('Output/fpkm.txt')

colnames(matrix) <- paste0('FPKM.', colnames(matrix))

write.table(matrix, 'Output/merged_result/genes_expression.txt', row.names = T, quote = F, sep = '\t')

#####
test <- read.delim('ref_transcripts_expression.txt')

matrix$gene_id <- rownames(matrix)

ggvenn(list(matrix=matrix$gene_id, test=test$gene_id))

final <- left_join(test[,1:10], matrix)

colnames(final)

cov <- final[,-c(1:10)]

colnames(cov) <- gsub('FPKM','cov',colnames(cov))

final2 <- cbind(final, cov)

write.table(final2, 'Output/merged_result/transcripts_expression.txt', row.names = T, quote = F, sep = '\t')





