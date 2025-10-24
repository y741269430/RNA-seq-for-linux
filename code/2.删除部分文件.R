# =============================================
# 删除diff_result里面关于 transcript 的文件 ####
# =============================================
setwd(r"{D:\R work\Output\}")

# 设置主文件夹路径
base_folder <- "COND1"

group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

# 定义要匹配的文件模式
patterns <- c(
  "*_edger_result_gene.DESeq2.R.log",  # 模式1
  "*_transcript*",                     # 模式2：任何包含 transcript 的文件
  "*_edger_result_gene.DESeq2.DE_results"
)

# 遍历每个子文件夹
for (subfolder in group) {
  # 构造 diff_result 文件夹的完整路径
  diff_result_path <- file.path(base_folder, subfolder, "diff_result")
  
  # 检查该 diff_result 文件夹是否存在
  if (!dir.exists(diff_result_path)) {
    cat("⚠️  diff_result 文件夹不存在: ", diff_result_path, "\n")
    next
  }
  
  cat("🔍 正在检查文件夹:", diff_result_path, "\n")
  
  # 遍历每个文件匹配模式
  for (pattern in patterns) {
    # 查找匹配的文件
    matched_files <- list.files(
      path = diff_result_path,
      pattern = pattern,
      full.names = TRUE
    )
    
    if (length(matched_files) == 0) {
      cat("   🚫 未找到匹配模式 '", pattern, "' 的文件\n", sep = "")
    } else {
      cat("   ✅ 找到", length(matched_files), "个匹配模式 '", pattern, "' 的文件\n", sep = "")
      
      # 遍历并删除每一个匹配的文件
      for (file_path in matched_files) {
        if (file.exists(file_path)) {
          file.remove(file_path)
          cat("     🗑️ 已删除文件:", file_path, "\n")
        } else {
          cat("     ⚠️ 文件不存在（可能已被删除）:", file_path, "\n")
        }
      }
    }
  }
}

# =============================================
# 删除mRNA_result里面匹配以下规则的文件 ####
# =============================================
setwd(r"{D:\R work\Output\}")

# 设置主文件夹路径
base_folder <- "COND1"

group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

# 定义要匹配的文件模式
patterns_to_delete <- c(
  "*_genes_annotation.txt",             # 类型 1
  "*_genes_annotation_volcano.log",     # 类型 2
  "*_transcripts_annotation.txt"        # 类型 3
)

# 遍历每个子文件夹
for (mRNA_folder in group) {
  # 构造 mRNA_result 文件夹的完整路径
  result_path <- file.path(base_folder, mRNA_folder, "mRNA_result")
  
  # 检查该 mRNA_result 文件夹是否存在
  if (!dir.exists(result_path)) {
    cat("⚠️  mRNA_result 文件夹不存在: ", result_path, "\n")
    next  # 跳过这个子文件夹
  }
  
  cat("🔍 正在检查文件夹:", result_path, "\n")
  
  # 遍历每一个要删除的文件模式
  for (pattern in patterns_to_delete) {
    # 查找所有匹配该模式的文件（返回完整路径）
    matched_files <- list.files(
      path = result_path,
      pattern = pattern,
      full.names = TRUE
    )
    
    if (length(matched_files) == 0) {
      cat("   🚫 未找到匹配模式 '", pattern, "' 的文件\n", sep = "")
    } else {
      cat("   ✅ 找到", length(matched_files), "个匹配模式 '", pattern, "' 的文件\n", sep = "")
      
      # 筛选：排除以 add_ 开头的文件
      files_to_delete <- matched_files[!grepl("^add_", basename(matched_files))]
      
      if (length(files_to_delete) == 0) {
        cat("   ℹ️  所有匹配文件均以 add_ 开头，没有需要删除的。\n")
      } else {
        cat("   🗑️  将删除", length(files_to_delete), "个文件（排除以 add_ 开头的）\n")
        
        # 逐个删除
        for (file_path in files_to_delete) {
          if (file.exists(file_path)) {
            file.remove(file_path)
            cat("     🗑️ 已删除文件:", file_path, "\n")
          } else {
            cat("     ⚠️ 文件不存在（可能已被删除）:", file_path, "\n")
          }
        }
      }
    }
  }
}


# =============================================
# 删除gsea_result下的部分文件 ####
# =============================================
setwd(r"{D:\R work\Output\}")

base_folder <- "COND1"
group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

# ---- (1) gsea_result/GO/ 文件夹中的文件 ----
go_fixed <- c("GSEA_GO.log", "GSEA_GO.sh")
go_wildcard <- c("*.Gsea.enrichment.GO.rds", "*.Gsea.enrichment.GO.txt")

# ---- (2) gsea_result/KEGG/ 文件夹中的文件 ----
kegg_fixed <- c("GSEA_KEGG.log", "GSEA_KEGG.sh")
kegg_wildcard <- c("*.Gsea.enrichment.KEGG.rds", "*.Gsea.enrichment.KEGG.txt")

# ---- (3) gsea_result/ 文件夹中的文件 ----
gsea_root_wildcard <- c("*.Gsea.enrichment.cls", "*.Gsea.enrichment.Expression.txt")

for (subfolder in group) {
  # ===== (1) 处理 gsea_result/GO/ 文件夹 =====
  go_path <- file.path(base_folder, subfolder, "gsea_result", "GO")
  if (dir.exists(go_path)) {
    cat("🔍 检查 GO 文件夹:", go_path, "\n")
    for (fname in go_fixed) {
      file_path <- file.path(go_path, fname)
      if (file.exists(file_path)) {
        file.remove(file_path)
        cat("     🗑️ 已删除文件:", file_path, "\n")
      } else {
        cat("     🚫 文件不存在:", file_path, "\n")
      }
    }
    for (pattern in go_wildcard) {
      matched <- list.files(path = go_path, pattern = pattern, full.names = TRUE)
      if (length(matched) > 0) {
        for (f in matched) {
          file.remove(f)
          cat("     🗑️ 已删除通配文件:", f, "\n")
        }
      } else {
        cat("     🚫 未找到匹配通配符 '", pattern, "' 的文件\n", sep = "")
      }
    }
  } else {
    cat("⚠️  GO 文件夹不存在:", go_path, "\n")
  }
  
  # ===== (2) 处理 gsea_result/KEGG/ 文件夹 =====
  kegg_path <- file.path(base_folder, subfolder, "gsea_result", "KEGG")
  if (dir.exists(kegg_path)) {
    cat("🔍 检查 KEGG 文件夹:", kegg_path, "\n")
    for (fname in kegg_fixed) {
      file_path <- file.path(kegg_path, fname)
      if (file.exists(file_path)) {
        file.remove(file_path)
        cat("     🗑️ 已删除文件:", file_path, "\n")
      } else {
        cat("     🚫 文件不存在:", file_path, "\n")
      }
    }
    for (pattern in kegg_wildcard) {
      matched <- list.files(path = kegg_path, pattern = pattern, full.names = TRUE)
      if (length(matched) > 0) {
        for (f in matched) {
          file.remove(f)
          cat("     🗑️ 已删除通配文件:", f, "\n")
        }
      } else {
        cat("     🚫 未找到匹配通配符 '", pattern, "' 的文件\n", sep = "")
      }
    }
  } else {
    cat("⚠️  KEGG 文件夹不存在:", kegg_path, "\n")
  }
  
  # ===== (3) 处理 gsea_result/ 文件夹 =====
  gsea_root_path <- file.path(base_folder, subfolder, "gsea_result")
  if (dir.exists(gsea_root_path)) {
    cat("🔍 检查 gsea_result 根文件夹:", gsea_root_path, "\n")
    for (pattern in gsea_root_wildcard) {
      matched <- list.files(path = gsea_root_path, pattern = pattern, full.names = TRUE)
      if (length(matched) > 0) {
        for (f in matched) {
          file.remove(f)
          cat("     🗑️ 已删除通配文件:", f, "\n")
        }
      } else {
        cat("     🚫 未找到匹配通配符 '", pattern, "' 的文件\n", sep = "")
      }
    }
  } else {
    cat("⚠️  gsea_result 根文件夹不存在:", gsea_root_path, "\n")
  }
}


# =============================================
# 删除gsea_result里面关于 my_analysis.Gsea 的文件 （权限不足，手动整理）####
# =============================================

setwd(r"{D:\R work\Output\}")  

base_folder <- "COND1"
group <- c('SUSVSCtrl', 'CtrlPVSCtrl', 'SUSPVSCtrlP', 'SUSPVSSUS')

# 定义要删除的文件名模式
pattern_to_delete <- "my_analysis.Gsea*"  # 匹配如 my_analysis.Gsea*

for (subfolder in group) {
  # 构造 gsea_result/GO/ 和 gsea_result/KEGG/ 的路径
  go_path <- file.path(base_folder, subfolder, "gsea_result", "GO")
  kegg_path <- file.path(base_folder, subfolder, "gsea_result", "KEGG")
  
  # --------------------------
  # 处理 GO 文件夹
  # --------------------------
  if (dir.exists(go_path)) {
    cat("🔍 检查 GO 文件夹:", go_path, "\n")
    matched_go_files <- list.files(
      path = go_path,
      pattern = pattern_to_delete,
      full.names = TRUE
    )
    
    if (length(matched_go_files) == 0) {
      cat("   🚫 GO 文件夹中没有找到匹配 '", pattern_to_delete, "' 的文件\n", sep = "")
    } else {
      cat("   ✅ GO 文件夹中找到", length(matched_go_files), "个匹配文件，准备删除...\n", sep = "")
      for (file_path in matched_go_files) {
        if (file.exists(file_path)) {
          file.remove(file_path)
          cat("     🗑️ 已删除文件:", file_path, "\n")
        } else {
          cat("     ⚠️ 文件不存在（可能已删除）:", file_path, "\n")
        }
      }
    }
  } else {
    cat("⚠️  GO 文件夹不存在:", go_path, "\n")
  }
  
  # --------------------------
  # 处理 KEGG 文件夹
  # --------------------------
  if (dir.exists(kegg_path)) {
    cat("🔍 检查 KEGG 文件夹:", kegg_path, "\n")
    matched_kegg_files <- list.files(
      path = kegg_path,
      pattern = pattern_to_delete,
      full.names = TRUE
    )
    
    if (length(matched_kegg_files) == 0) {
      cat("   🚫 KEGG 文件夹中没有找到匹配 '", pattern_to_delete, "' 的文件\n", sep = "")
    } else {
      cat("   ✅ KEGG 文件夹中找到", length(matched_kegg_files), "个匹配文件，准备删除...\n", sep = "")
      for (file_path in matched_kegg_files) {
        if (file.exists(file_path)) {
          file.remove(file_path)
          cat("     🗑️ 已删除文件:", file_path, "\n")
        } else {
          cat("     ⚠️ 文件不存在（可能已删除）:", file_path, "\n")
        }
      }
    }
  } else {
    cat("⚠️  KEGG 文件夹不存在:", kegg_path, "\n")
  }
}



# =============================================
# 删除.rds文件
# =============================================
setwd(r"{D:\R work\Output1\}")

# 设置主文件夹名称
base_folder <- "COND1"

# 构造 COND1 文件夹的完整路径
cond1_path <- file.path(base_folder)

# 递归查找所有 .rds 文件（包括所有子文件夹中的 .rds）
rds_files <- list.files(
  path = cond1_path,         # 检索 COND1 文件夹
  pattern = "\\.rds$",       # 匹配以 .rds 结尾的文件
  full.names = TRUE,         # 返回完整路径
  recursive = TRUE           # 递归查找子文件夹
)

# 判断是否找到 .rds 文件
if (length(rds_files) == 0) {
  cat("   🚫 未找到任何 .rds 文件\n")
} else {
  cat("   ✅ 找到", length(rds_files), "个 .rds 文件，准备删除...\n", sep = "")
  
  # 7. 遍历每一个 .rds 文件，删除它
  for (rds_file in rds_files) {
    if (file.exists(rds_file)) {
      # 可选：打印将要删除的文件路径
      cat("     🗑️ 删除文件:", rds_file, "\n")
      
      # 删除文件
      file.remove(rds_file)
      
      # 可选：再次检查是否真的删除了（调试用）
      if (!file.exists(rds_file)) {
        # cat("     ✅ 已成功删除\n")
      } else {
        cat("     ⚠️ 删除失败（文件可能被占用或权限不足）:", rds_file, "\n")
      }
    } else {
      cat("     ⚠️ 文件不存在（可能已被删除）:", rds_file, "\n")
    }
  }
}






