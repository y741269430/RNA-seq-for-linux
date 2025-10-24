
setwd(r"{D:\R work\Output1\}")  

base_folder <- "COND1"
group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

data <- read.delim('COND1/*/mRNA_result/add_*_genes_annotation.txt')
write.xlsx(data, 'COND1/*/mRNA_result/DEG_*_genes_annotation.xlsx')

能否批量对这些文件进行txt转xlsx
COND1\SUSPVSSUS\mRNA_result

#### 批量对mRNA_result里面的txt，以及func_result里面的txt，转成xlsx


# =============================================
# 1. 设置工作目录和主文件夹
# =============================================
setwd(r"{D:\R work\Output1\}")

base_folder <- "COND1"
group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

# =============================================
# 2. 遍历每个 group（子文件夹）
# =============================================
for (subfolder in group) {
  # 构造 mRNA_result 文件夹路径
  txt_folder_path <- file.path(base_folder, subfolder, "mRNA_result")
  
  # 检查该文件夹是否存在
  if (!dir.exists(txt_folder_path)) {
    cat("⚠️ mRNA_result 文件夹不存在:", txt_folder_path, "\n")
    next
  }
  
  cat("🔍 正在处理子文件夹:", subfolder, "路径:", txt_folder_path, "\n")
  
  # 查找所有符合 add_*_genes_annotation.txt 的文件
  txt_files <- list.files(
    path = txt_folder_path,
    pattern = "^add_.*_genes_annotation\\.txt$",  
    full.names = TRUE
  )
  
  if (length(txt_files) == 0) {
    cat("   🚫 未找到匹配 add_*_genes_annotation.txt 的文件\n")
    next
  }
  
  cat("   ✅ 找到", length(txt_files), "个匹配的 .txt 文件\n", sep = "")
  
  # 遍历每个 txt 文件，读取并转为 xlsx
  for (txt_file in txt_files) {
    # 读取 txt 文件（假设是制表符分隔）
    data <- tryCatch({
      read.delim(txt_file, stringsAsFactors = FALSE)
    }, error = function(e) {
      cat("     ❌ 读取文件失败:", txt_file, "-", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(data)) next  # 读取失败则跳过

    xlsx_file <- sub('add', 'DEG', sub("\\.txt$", ".xlsx", txt_file))

    # 写入 xlsx 文件
    tryCatch({
      write.xlsx(data, xlsx_file)
      cat("     🔄 已转换: ", txt_file, " → ", xlsx_file, "\n", sep = "")
    }, error = function(e) {
      cat("     ❌ 写入 xlsx 失败:", xlsx_file, "-", e$message, "\n")
    })
  }
}


# =============================================
# 删除原来的 add_*_genes_annotation.txt 文件
# =============================================

# 1. 遍历每个 group（子文件夹）
for (subfolder in group) {
  txt_folder_path <- file.path(base_folder, subfolder, "mRNA_result")
  
  if (!dir.exists(txt_folder_path)) {
    cat("⚠️ mRNA_result 文件夹不存在:", txt_folder_path, "\n")
    next
  }
  
  cat("🔍 正在处理子文件夹:", subfolder, "路径:", txt_folder_path, "\n")
  
  # 2. 查找所有 add_*_genes_annotation.txt 文件
  txt_files <- list.files(
    path = txt_folder_path,
    pattern = "^add_.*_genes_annotation\\.txt$",
    full.names = TRUE
  )
  
  if (length(txt_files) == 0) {
    cat("   🚫 未找到匹配的 .txt 文件\n")
    next
  }
  
  cat("   ✅ 找到", length(txt_files), "个 .txt 文件，准备删除\n", sep = "")
  
  # 3. 遍历每个 txt 文件，尝试删除
  for (txt_file in txt_files) {
    # 4. 检查文件是否存在
    if (file.exists(txt_file)) {
      # 5. 删除原 txt 文件
      file.remove(txt_file)
      cat("     🗑️ 已删除原文件:", txt_file, "\n")
    } else {
      cat("     ⚠️ 原文件不存在，无法删除:", txt_file, "\n")
    }
  }
}




# =============================================
# 1. 设置工作目录和主文件夹
# =============================================
setwd(r"{D:\R work\Output\}") 

base_folder <- "COND1"
group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

# =============================================
# 2. 定义三类要转换的文件模式（在 func_result/ 目录下）
# =============================================

# 每个 pattern 是在 COND1/<subfolder>/func_result/ 下的文件名匹配规则
patterns <- list(
  GO_Enrichment = "COND1.*.GO_Enrichment.txt",
  GO_Enrichment_gene = "COND1.*.GO_Enrichment_gene.txt",
  KEGG_Enrichment_gene = "COND1.*.KEGG_Enrichment_gene.txt"
)

# =============================================
# 3. 遍历每个 group（子文件夹），查找并转换文件
# =============================================
for (subfolder in group) {
  cat("🔍 正在处理子文件夹:", subfolder, "\n")
  
  func_result_path <- file.path("COND1", subfolder, "func_result")
  
  if (!dir.exists(func_result_path)) {
    cat("   ⚠️ func_result 文件夹不存在:", func_result_path, "\n")
    next
  }
  
  # 遍历每一类文件模式
  for (type in names(patterns)) {
    pattern <- patterns[[type]]
    # 查找匹配的 .txt 文件
    txt_files <- list.files(
      path = func_result_path,
      pattern = pattern,
      full.names = TRUE
    )
    
    if (length(txt_files) == 0) {
      cat("   🚫 未找到匹配模式 '", pattern, "' 的文件\n", sep = "")
      next
    }
    
    cat("   ✅ 找到", length(txt_files), "个 '", type, "' 文件\n", sep = "")
    
    for (txt_file in txt_files) {
      # 构造对应的 xlsx 文件名（将 .txt 替换为 .xlsx）
      xlsx_file <- sub("\\.txt$", ".xlsx", txt_file)
      
      # 读取 txt 文件（假设是制表符分隔，可根据实际调整 sep = "," 等）
      data <- tryCatch({
        read.delim(txt_file, stringsAsFactors = FALSE)
      }, error = function(e) {
        cat("     ❌ 读取文件失败:", txt_file, "-", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(data)) next  # 读取失败则跳过
      
      # 写入 xlsx 文件
      tryCatch({
        write.xlsx(data, xlsx_file)
        cat("     🔄 已转换: ", txt_file, " → ", xlsx_file, "\n", sep = "")
      }, error = function(e) {
        cat("     ❌ 写入 xlsx 失败:", xlsx_file, "-", e$message, "\n")
      })
    }
  }
}
