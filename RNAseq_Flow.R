library(clusterProfiler)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(org.Mm.eg.db)
library(pheatmap)
library(openxlsx)
library(enrichplot)
library(factoextra)
library(ggthemes)
library(cowplot)

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
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)  
  nor <- DESeq2::counts(dds, normalized = T)
  nor2 <- as.data.frame(nor)
}
myPCA <- function(x, grouplist, label_option = "ind", title = NULL) {
  factoextra::fviz_pca_ind(prcomp(t(x),scale. = T), 
                           repel = TRUE, pointsize = 5, 
                           palette = "jco", label = ifelse(label_option == "none", "none", "ind"), mean.point = F, 
                           col.ind = grouplist) + 
    ggplot2::coord_fixed(1) + 
    ggplot2::labs(title = element_blank())+ 
    ggplot2::theme_bw()+
    ggplot2::theme(text = element_text(family = "sans"),
                   plot.title = element_text(hjust = 0.5),
                   axis.title = element_text(size = 15),
                   axis.text = element_text(size = 15),
                   legend.title = element_blank(), 
                   legend.text = element_text(face = "bold", size = 15),
                   legend.position = "bottom") + 
    xlim(-250, 250) + ylim(-250, 250) + ggtitle(title)
}
myDESeq2 <- function(x, ct, tr, id = 'gene_id'){
  # raw counts
  localgeneid <- read.delim('rmdup_20251102.txt')
  
  condition <- factor(c(rep("ctrl", ct), rep("treatment", tr)), 
                      levels = c("ctrl", "treatment"))
  
  colData <- data.frame(row.names = colnames(x), condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                        colData = colData,
                                        design = ~condition)
  # keep <- rowSums(DESeq2::counts(dds)) >= 100
  # dds <- dds[keep,]
  
  dds <- DESeq2::DESeq(dds)    
  nor <- DESeq2::counts(dds, normalized = T)
  
  DEG_DESeq2 <- merge(data.frame(results(dds, contrast = c("condition", "treatment", "ctrl"))), nor, by = "row.names")
  DEG_DESeq2 <- na.omit(DEG_DESeq2)
  colnames(DEG_DESeq2)[1] <- id
  colnames(DEG_DESeq2)[which(colnames(DEG_DESeq2) == 'pvalue')] = 'pvalue'
  colnames(DEG_DESeq2)[which(colnames(DEG_DESeq2) == 'log2FoldChange')] = 'log2FC'
  DEG <- merge(localgeneid, DEG_DESeq2, id)
  
  return(DEG)
}
myVol <- function(expr, lg = log2(1.5), title = NULL, sort_by_fc = TRUE){
  expr$log10pvalue <- -log10(expr$pvalue)
  expr <- expr[!is.infinite(expr$log10pvalue),]
  expr$Group <- "None-significant"
  expr$Group[which((expr$pvalue < 0.05) & (expr$log2FC > lg))] <- "Up-regulated"
  expr$Group[which((expr$pvalue < 0.05) & (expr$log2FC < -lg))] <- "Down-regulated"
  table(expr$Group)
  
  expr$label = ""
  
  # 根据参数选择排序方式
  if (sort_by_fc) {
    expr <- expr[order(expr$log2FC, decreasing = TRUE), ]
    upgenes <- head(expr$SYMBOL[which(expr$Group == "Up-regulated")], 10)
    downgenes <- tail(expr$SYMBOL[which(expr$Group == "Down-regulated")], 10)
    
  } else {
    expr <- expr[order(expr$pvalue), ]
    upgenes <- head(expr$SYMBOL[which(expr$Group == "Up-regulated")], 10)
    downgenes <- head(expr$SYMBOL[which(expr$Group == "Down-regulated")], 10)
  }
  
  top10deg <- c(as.character(upgenes), as.character(downgenes))
  expr$label[match(top10deg, expr$SYMBOL)] <- top10deg
  
  expr <- subset(expr, log10pvalue < quantile(log10pvalue, 0.9999))
  
  # 计算对称范围的最大值
  symmetric_range <- max(abs(min(expr$log2FC)), abs(max(expr$log2FC)))+1
  
  # 对称范围的最小值和最大值
  min_symmetric <- -symmetric_range
  max_symmetric <- symmetric_range
  
  p <- ggpubr::ggscatter(expr, 
                         x = "log2FC", 
                         y = "log10pvalue", 
                         color = "Group", 
                         palette = c("#354E6B", "grey", "#A64036"), 
                         size = 1,
                         label = expr$label,
                         font.label = c(14, "bold"), 
                         font.family = "sans",
                         repel = T, 
                         xlab = "log2FoldChange", 
                         ylab = "-log10(pvalue)") + 
    theme_base() + 
    geom_hline(yintercept = 1.3, linetype = "dashed") + 
    geom_vline(xintercept = c(-lg,lg), linetype = "dashed") +
    theme(legend.title = element_blank(),
          legend.position = "top",
          text = element_text(family = "sans")) + 
    ggtitle(title) + 
    xlim(min_symmetric, max_symmetric)

  p + guides(color = guide_legend(override.aes = list(label = "")))
} 
myHeat <- function(x, grouplist, show_rownames = F, 
                   angle_col = 45, cluster_cols = F, cluster_rows = T, main = NA){
  # 输入一个标准化的矩阵，进行热图绘制，以grouplist为准，默认关闭行名，列名45度
  # 输出是一个热图
  # 保存热图一定要加上pheatmap的包名: pheatmap::
  
  group <- data.frame(grouplist)
  rownames(group) <- colnames(x)
  plot <- pheatmap::pheatmap(log10(x + 1), 
                             scale = "row", 
                             clustering_distance_rows = "correlation",
                             #color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(256)),
                             color = colorRampPalette(c("navy", "white", "firebrick3"))(256),
                             # color = colorRampPalette(c('#03045E', '#0077B6', '#00B4D8', '#90E0EF', '#CAF0F8',
                             #                            '#FAE0E4', '#F9BEC7', '#FF99AC', '#FF7096', '#FF477E'))(256),
                             fontsize = 6, 
                             fontfamily="sans",
                             display_numbers = F,
                             border_color = NA,
                             gaps_row = F,
                             cluster_cols = cluster_cols, 
                             cluster_rows = cluster_rows,
                             treeheight_row = 50, treeheight_col = 50,
                             cellwidth = NA, 
                             cellheight = NA,
                             legend = T,
                             show_rownames = show_rownames, 
                             show_colnames = T,
                             annotation_legend = T,
                             annotation_names_col = F,
                             annotation_col = group, 
                             angle_col = angle_col, 
                             main = main,
                             silent = T) 
  plot_grid(plot$gtable)
}
my_enrichKEGG <- function(gene, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH", universe, minGSSize = 10, maxGSSize = 500, 
                          qvalueCutoff = 0.2, use_internal_data = FALSE) {
  species <- clusterProfiler:::organismMapper(organism)
  if (use_internal_data) {KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db(species)}
  else {
    astk_dir <- file.path("./")
    kegg_data_name <- paste(organism, lubridate::today(),"kegg","RData", sep = ".")
    kegg_data_path <- file.path(astk_dir, kegg_data_name)
    if (file.exists(kegg_data_path)){
      print("load local catch")
      load(kegg_data_path)
    }else {
      print("online downloading... ")
      KEGG_DATA <- clusterProfiler:::prepare_KEGG(organism,"KEGG", "kegg")
      save(KEGG_DATA, file = kegg_data_path) 
    }; print(KEGG_DATA)
    
  }
  
  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff, 
                                             pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
                                             maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = KEGG_DATA)
  if (is.null(res)) 
    return(res)
  res@ontology <- "KEGG"
  res@organism <- species
  res@keytype <- keyType
  return(res)
}
localgeneid <- read.delim('rmdup_20251102.txt')
