# =============================================
# 横向排列的分面火山图 
# 日期: 2025-08-14
# 用户: buyu-ahau
# =============================================

# 加载必要的包
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(cowplot)
library(tidyr)

# 设置工作目录
setwd("/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/scRNA-seq数据")

# 定义文件路径
de_results_directory <- "差异表达分析/基因符号_DE"
output_directory <- "差异表达分析/火山图_灰色小FC"

# 创建输出目录
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# =============================================
# 步骤1：收集和准备数据
# =============================================

cat("收集所有细胞类型的数据...\n")

# 获取所有细胞类型
cell_types <- list.dirs(de_results_directory, full.names = FALSE, recursive = FALSE)
cell_types <- cell_types[!grepl("统计", cell_types)]  # 排除可能的非细胞类型目录

# 存储所有数据
all_data <- data.frame()
extreme_counts <- c()

# 处理每种细胞类型
for(cell_type in cell_types) {
  cat("处理", cell_type, "的数据...\n")
  
  # 读取显著差异基因数据
  sig_file <- file.path(de_results_directory, cell_type, 
                      paste0(cell_type, "_significant_DE_genes_with_symbols.csv"))
  
  if(file.exists(sig_file)) {
    data <- fread(sig_file)
    
    # 分离极显著基因
    extreme_genes <- data %>% filter(p_val_adj <= 1e-300)
    regular_genes <- data %>% filter(p_val_adj > 1e-300)
    
    # 记录极显著基因数量
    extreme_counts[cell_type] <- nrow(extreme_genes)
    
    # 为regular_genes添加细胞类型信息
    regular_genes$cell_type <- cell_type
    
    # 添加到总数据集
    all_data <- bind_rows(all_data, regular_genes)
    
    # 输出统计
    cat("  总基因数:", nrow(data), "\n")
    cat("  常规p值基因数 (>1e-300):", nrow(regular_genes), "\n")
    cat("  极显著基因数 (<=1e-300):", nrow(extreme_genes), "\n\n")
  } else {
    cat("  未找到文件:", sig_file, "\n\n")
  }
}

# =============================================
# 步骤2：准备绘图数据
# =============================================

cat("准备绘图数据...\n")

# 计算-log10(p_val_adj)
all_data$neg_log10_padj <- -log10(all_data$p_val_adj)

# 定义新的调控状态分类 - 修改为只有高变化的基因有颜色
all_data$regulation <- case_when(
  # 显著上调且FC>1的基因
  all_data$p_val_adj < 0.05 & all_data$avg_log2FC >= 1 ~ "Highly Upregulated",
  
  # 显著下调且FC<-1的基因
  all_data$p_val_adj < 0.05 & all_data$avg_log2FC <= -1 ~ "Highly Downregulated",
  
  # 其余所有基因设为灰色（包括显著差异但FC小的基因和非显著基因）
  TRUE ~ "Not Significant/Small FC"
)

# 转换为因子并设置顺序
all_data$regulation <- factor(all_data$regulation, 
                            levels = c("Highly Upregulated", "Highly Downregulated", 
                                      "Not Significant/Small FC"))

# 转换细胞类型名称为展示名称
all_data <- all_data %>%
  mutate(
    display_name = case_when(
      cell_type == "B cells" ~ "B_cell",
      cell_type == "CD4 cells" ~ "CD4_cell",
      cell_type == "CD8 cells" ~ "CD8_cell",
      cell_type == "DC cells" ~ "DC_cell",
      cell_type == "Endothelial cells" ~ "Endothelial_cell",
      cell_type == "Epithelial cells" ~ "Epithelial_cell",
      cell_type == "Fibroblasts" ~ "Fibroblast",
      cell_type == "Macrophages" ~ "Macrophage",
      cell_type == "Neutrophils" ~ "Neutrophil",
      cell_type == "NK cells" ~ "NK-T_cell",
      cell_type == "Plasma cells" ~ "Plasma_cell",
      TRUE ~ as.character(cell_type)
    )
  )

# 确保所有11个细胞类型都被包括在自定义顺序中
custom_display_order <- c(
  "B_cell", "CD4_cell", "CD8_cell", "DC_cell",
  "Endothelial_cell", "Epithelial_cell", "Fibroblast",
  "Macrophage", "Neutrophil", "NK-T_cell", "Plasma_cell"
)

# 只保留实际存在的细胞类型
custom_display_order <- custom_display_order[custom_display_order %in% unique(all_data$display_name)]

# 设置细胞类型因子顺序
all_data$display_name <- factor(all_data$display_name, levels = custom_display_order)

# =============================================
# 步骤3：设置颜色方案
# =============================================

# 为每个细胞类型设置不同的颜色
cell_colors <- c(
  "B_cell" = "firebrick",
  "CD4_cell" = "skyblue3",
  "CD8_cell" = "forestgreen",
  "DC_cell" = "deepskyblue3",
  "Endothelial_cell" = "chartreuse4",
  "Epithelial_cell" = "chartreuse3",
  "Fibroblast" = "hotpink",
  "Macrophage" = "gold2",
  "Neutrophil" = "turquoise3",
  "NK-T_cell" = "mediumpurple",
  "Plasma_cell" = "steelblue"
)

# 为调控状态设置颜色 - 修改颜色方案
regulation_colors <- c(
  "Highly Upregulated" = "#E41A1C",       # 红色
  "Highly Downregulated" = "#377EB8",     # 蓝色
  "Not Significant/Small FC" = "#CCCCCC"  # 灰色
)

# =============================================
# 步骤4：创建单排显示的火山图
# =============================================

cat("\n创建横向排列的火山图...\n")

# 创建单行排列的火山图
single_row_plot <- ggplot(all_data, aes(x = neg_log10_padj, y = avg_log2FC)) +
  # 添加点
  geom_point(aes(color = regulation), size = 0.7, alpha = 0.8) +
  # 设置颜色
  scale_color_manual(values = regulation_colors, 
                     name = "Gene Regulation",
                     labels = c("Upregulated (FC ≥ 1)", 
                               "Downregulated (FC ≤ -1)", 
                               "Other (p>0.05 or |FC|<1)")) +
  # 水平参考线
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.5) +
  # 垂直参考线
  geom_vline(xintercept = -log10(0.05), 
             linetype = "dashed", color = "grey50", alpha = 0.5) +
  # 分面
  facet_wrap(~ display_name, nrow = 1, scales = "free_x") +
  # 标题和坐标轴标签
  labs(
    title = "Faceted Volcano Plot of DEGs (Day28 vs Day0)",
    subtitle = "Only genes with |FC| >= 1 are highlighted",
    x = "-log10 (Adjusted P-value)",
    y = "Average Log2 Fold Change"
  ) +
  # 固定Y轴范围
  coord_cartesian(ylim = c(-10, 10)) +
  # 美化主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey50", fill = NA),
    panel.spacing = unit(0.2, "lines"),
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 保存单行排列的图表
ggsave(file.path(output_directory, "Grey_Small_FC_Volcano_Row.png"), 
      plot = single_row_plot, width = 28, height = 8, dpi = 300, bg = "white")
ggsave(file.path(output_directory, "Grey_Small_FC_Volcano_Row.pdf"), 
      plot = single_row_plot, width = 28, height = 8, bg = "white")

# =============================================
# 步骤5：创建3行4列的火山图布局
# =============================================

cat("\n创建3行布局的火山图...\n")

# 每行显示的细胞类型数
cells_per_row <- 4

# 计算需要的行数
total_cells <- length(custom_display_order)
num_rows <- ceiling(total_cells / cells_per_row)

# 分配细胞类型到不同行
row_cell_types <- list()
for(i in 1:num_rows) {
  start_idx <- (i-1) * cells_per_row + 1
  end_idx <- min(i * cells_per_row, total_cells)
  if(start_idx <= total_cells) {
    row_cell_types[[i]] <- custom_display_order[start_idx:end_idx]
  }
}

# 为每行创建单独的图表
row_plots <- list()

for(i in 1:length(row_cell_types)) {
  # 获取当前行的细胞类型
  current_row_cells <- row_cell_types[[i]]
  
  # 筛选当前行的数据
  row_data <- all_data %>% filter(display_name %in% current_row_cells) %>%
    mutate(display_name = factor(display_name, levels = current_row_cells))
  
  # 创建当前行的火山图
  p <- ggplot(row_data, aes(x = neg_log10_padj, y = avg_log2FC)) +
    geom_point(aes(color = regulation), size = 0.7, alpha = 0.8) +
    scale_color_manual(values = regulation_colors, name = "Gene Regulation") +
    geom_hline(yintercept = -1, linetype = "dashed", color = "grey50", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50", alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50", alpha = 0.5) +
    facet_wrap(~ display_name, nrow = 1, scales = "free_x") +
    labs(
      x = if(i == length(row_cell_types)) "-log10 (Adjusted P-value)" else "",
      y = "Average Log2 Fold Change"
    ) +
    coord_cartesian(ylim = c(-10, 10)) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey50", fill = NA),
      panel.spacing = unit(0.2, "lines"),
      legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 8)
    )
  
  row_plots[[i]] <- p
}

# 添加标题
title <- ggdraw() + 
  draw_label("Faceted Volcano Plot of DEGs (Day28 vs Day0)", 
            fontface = "bold", size = 16, x = 0.5) +
  draw_label("Only genes with |FC| >= 1 are highlighted", 
            size = 12, x = 0.5, y = 0.3)

# 创建图例
legend_df <- data.frame(
  x = 1:3,
  y = rep(1, 3),
  regulation = factor(1:3, levels = 1:3, labels = levels(all_data$regulation))
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = regulation)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = regulation_colors,
    name = "Gene Regulation",
    labels = c("Upregulated (FC ≥ 1)", 
               "Downregulated (FC ≤ -1)", 
               "Other (p>0.05 or |FC|<1)")
  ) +
  theme_void() +
  theme(legend.position = "bottom")

# 提取图例
legend <- cowplot::get_legend(legend_plot)

# 组合标题、图表和图例
plot_list <- c(list(title), row_plots, list(legend))
heights <- c(0.1, rep(1, length(row_plots)), 0.1)

layout_plot <- plot_grid(
  plotlist = plot_list,
  ncol = 1,
  rel_heights = heights
)

# 保存3行布局的图表
ggsave(file.path(output_directory, "Grey_Small_FC_Volcano_Grid.png"), 
      plot = layout_plot, width = 20, height = 5 * num_rows + 2, dpi = 300, bg = "white")
ggsave(file.path(output_directory, "Grey_Small_FC_Volcano_Grid.pdf"), 
      plot = layout_plot, width = 20, height = 5 * num_rows + 2, bg = "white")

cat("\n处理完成！所有结果已保存到:", output_directory, "\n")