# 加载必要包
library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# --- 1. 参数设置 ---
# 直接使用聚类名称作为细胞类型
cell_types <- c("B_cells", "CD4_CD8_T_cells", "Endothelial_cells", "Epithelial_cells", 
                "Fibroblasts", "Macrophages", "Neutrophils", "NK_cells", "Plasma_cells")

# --- 路径设置 ---
base_path <- "/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq"
data_dir <- file.path(base_path, "2025.8.14")
output_dir <- file.path(base_path, "differential_analysis_results")
plot_dir <- file.path(output_dir, "violin_plots")
intermediate_dir <- file.path(output_dir, "intermediate_da_peaks") 
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2. 读取Seurat对象 ---
rda_path <- file.path(data_dir, "daa52225-obj.Rda")
cat("Loading Seurat object from:", rda_path, "\n"); flush.console()
load(rda_path) # 这将加载名为"obj"的对象
cat("Seurat object loaded successfully. Object size:", format(object.size(obj), units = "Gb"), "\n"); flush.console()

# --- 3. 检查并设置细胞类型标识 ---
# 检查元数据中是否已经有细胞类型注释
if(!"cell_type" %in% colnames(obj@meta.data)) {
  cat("WARNING: No 'cell_type' column found in metadata.\n")
  cat("Using seurat_clusters directly for analysis.\n")
  
  # 直接使用seurat_clusters作为cell_type
  obj$cell_type <- obj$seurat_clusters
  
  cat("Created cell_type column from seurat_clusters.\n")
  cat("Cell types:", paste(sort(unique(obj$cell_type)), collapse=", "), "\n")
}

# --- 4. 设置Assay和分群 ---
assay_name <- if ("ATAC" %in% Assays(obj)) "ATAC" else "peaks"
DefaultAssay(obj) <- assay_name
cat("Set default assay to:", assay_name, "\n")

# 设置计数变量
latent_var <- if ("nCount_ATAC" %in% colnames(obj@meta.data)) "nCount_ATAC" else "nCount_peaks"
cat("Using", latent_var, "as latent variable for regression\n")

# 设置分群标识为cell_type
Idents(obj) <- obj$cell_type
cat("Set identities to cell_type. Available levels:", paste(levels(Idents(obj)), collapse=", "), "\n")

# --- 5. 循环进行差异可及性分析 ---
cat("Starting pairwise differential accessibility analysis...\n"); flush.console()

# 获取实际可用的细胞类型
available_cell_types <- levels(Idents(obj))
cat("Available cell types for analysis:", paste(available_cell_types, collapse=", "), "\n")

# 只使用可用的细胞类型进行比较
for (i in 1:(length(available_cell_types) - 1)) {
  for (j in (i + 1):length(available_cell_types)) {
    group1 <- available_cell_types[i]
    group2 <- available_cell_types[j]
    
    # 断点续跑功能: 检查中间文件是否已存在
    intermediate_file_path <- file.path(intermediate_dir, paste0("DApeaks_raw_", group1, "_vs_", group2, ".rds"))
    if (file.exists(intermediate_file_path)) {
        cat("Result for", group1, "vs", group2, "already exists. Skipping.\n"); flush.console()
        next # 如果文件存在，直接跳到下一次循环
    }
    
    # 检查每个组中是否有足够的细胞
    cells_group1 <- WhichCells(obj, idents = group1)
    cells_group2 <- WhichCells(obj, idents = group2)
    
    if (length(cells_group1) < 10 || length(cells_group2) < 10) {
      cat("Skipping comparison:", group1, "vs", group2, "- not enough cells in one or both groups.\n"); flush.console()
      cat("Cells in", group1, ":", length(cells_group1), ", Cells in", group2, ":", length(cells_group2), "\n"); flush.console()
      next
    }
    
    cat("Comparing:", group1, "vs", group2, "\n"); flush.console()
    
    # 内存优化: 在循环内部预先筛选数据
    cat("Subsetting Seurat object for comparison...\n"); flush.console()
    cells_to_keep <- c(cells_group1, cells_group2)
    
    # 使用tryCatch捕获可能的错误
    tryCatch({
      # 创建临时子集用于分析
      temp_seurat_subset <- subset(obj, cells = cells_to_keep)
      
      # 在子集上运行 FindMarkers
      da_peaks <- FindMarkers(
        object = temp_seurat_subset,
        ident.1 = group1,
        ident.2 = group2,
        test.use = 'LR',
        latent.vars = latent_var,
        min.pct = 0.05
      )
      
      # 保存原始的、未经筛选的 FindMarkers 中间结果
      saveRDS(da_peaks, file = intermediate_file_path)
      cat("Saved raw intermediate results to:", intermediate_file_path, "\n"); flush.console()
      
      # 筛选、注释并保存最终结果
      sig_peaks <- da_peaks[da_peaks$p_val_adj < 0.05, ]
      
      if (nrow(sig_peaks) == 0) {
          cat("No significant DA peaks found for", group1, "vs", group2, ". Skipping annotation.\n"); flush.console()
          rm(temp_seurat_subset, da_peaks)
          gc() 
          next
      }
      
      cat("Found", nrow(sig_peaks), "significant differential peaks.\n"); flush.console()
      cat("Annotating significant peaks...\n"); flush.console()
      sig_peaks$peak <- rownames(sig_peaks)
      
      # 使用ClosestFeature进行注释（已经测试过，这个函数可以正常工作）
      closest_features_all <- ClosestFeature(obj, regions = rownames(sig_peaks))
      annotated_sig_peaks <- cbind(sig_peaks, closest_features_all[, c("gene_name", "closest_region", "distance")])
      
      outcsv_all <- file.path(output_dir, paste0("DApeaks_", group1, "_vs_", group2, "_allsig_annotated.csv"))
      write.csv(annotated_sig_peaks, outcsv_all, row.names = FALSE)
      cat("Saved all", nrow(annotated_sig_peaks), "significant peaks to:", outcsv_all, "\n"); flush.console()
      
      # 优化Top peaks选择和合并
      n_pos <- min(50, nrow(annotated_sig_peaks[annotated_sig_peaks$avg_log2FC > 0, ]))
      n_neg <- min(50, nrow(annotated_sig_peaks[annotated_sig_peaks$avg_log2FC < 0, ]))
      
      top_50_pos <- head(annotated_sig_peaks[order(-annotated_sig_peaks$avg_log2FC), ], n_pos)
      top_50_neg <- head(annotated_sig_peaks[order(annotated_sig_peaks$avg_log2FC), ], n_neg)
      top_50_df <- rbind(top_50_pos, top_50_neg)
      
      if (nrow(top_50_df) > 0) {
          outcsv_top50 <- file.path(output_dir, paste0("DApeaks_", group1, "_vs_", group2, "_top50_annotated.csv"))
          write.csv(top_50_df, outcsv_top50, row.names = FALSE)
          cat("Saved top", nrow(top_50_df), "peaks to:", outcsv_top50, "\n"); flush.console()
      }
      
      # 批量绘制小提琴图
      cat("Generating violin plots for top peaks...\n"); flush.console()
      n_pos_plot <- min(10, nrow(top_50_df[top_50_df$avg_log2FC > 0, ]))
      n_neg_plot <- min(10, nrow(top_50_df[top_50_df$avg_log2FC < 0, ]))
      
      plot_peaks_pos <- head(top_50_df[top_50_df$avg_log2FC > 0, ], n_pos_plot)
      plot_peaks_neg <- head(top_50_df[top_50_df$avg_log2FC < 0, ], n_neg_plot)
      plot_peaks_df <- rbind(plot_peaks_pos, plot_peaks_neg)
      
      if (nrow(plot_peaks_df) > 0) {
        # 使用tryCatch捕获绘图错误
        tryCatch({
          plot_peaks <- plot_peaks_df$peak
          vln_plot <- VlnPlot(
            obj,
            features = plot_peaks,
            idents = c(group1, group2), # 只绘制比较的两个组
            group.by = "cell_type",
            pt.size = 0,
            ncol = 2
          ) + labs(title = paste(group1, "vs", group2, "(Top", nrow(plot_peaks_df), "Peaks)"))
          pdf_path <- file.path(plot_dir, paste0("Vln_", group1, "_vs_", group2, "_top", nrow(plot_peaks_df), "_DApeaks.pdf"))
          ggsave(pdf_path, vln_plot, width = 12, height = 5 * ceiling(length(plot_peaks) / 2), limitsize = FALSE)
          cat("Saved violin plot to:", pdf_path, "\n"); flush.console()
        }, error = function(e) {
          cat("WARNING: Failed to generate violin plot:", conditionMessage(e), "\n")
          cat("Continuing without plotting.\n"); flush.console()
        })
      } else {
        cat("No peaks available for plotting.\n"); flush.console()
      }
      
    }, error = function(e) {
      cat("ERROR in processing", group1, "vs", group2, ":", conditionMessage(e), "\n")
      cat("Skipping this comparison and moving to the next.\n"); flush.console()
    })
    
    # 内存优化: 在每次循环结束时，清理临时对象
    if (exists("temp_seurat_subset")) rm(temp_seurat_subset)
    if (exists("da_peaks")) rm(da_peaks)
    if (exists("sig_peaks")) rm(sig_peaks)
    if (exists("annotated_sig_peaks")) rm(annotated_sig_peaks)
    if (exists("top_50_df")) rm(top_50_df)
    if (exists("plot_peaks_df")) rm(plot_peaks_df)
    if (exists("vln_plot")) rm(vln_plot)
    gc()
    cat("Finished comparison and cleaned memory.\n\n"); flush.console()
  }
}

# 保存分析完成的时间和结果路径
end_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cat("分析完成时间:", end_time, "\n")
cat("全部分析完成！结果保存在:", output_dir, "\n"); flush.console()