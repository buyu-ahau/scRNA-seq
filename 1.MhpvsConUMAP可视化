# =============================================
# 感染组vs对照组UMAP可视化 - 自定义颜色
# =============================================

# 加载必要的库
library(Seurat)
library(ggplot2)
library(dplyr)
library(qs)
library(patchwork)

# 设置工作目录并加载数据
setwd("/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/signac/流程1")
rna_data <- qs::qread("rna_final_annotated_new.qs")

# 创建输出目录
dir.create("感染vs对照组比较", showWarnings = FALSE)

# =============================================
# 准备对照组和感染组数据
# =============================================

# 从细胞名称创建组别变量
print("从细胞名称创建组别变量...")
cell_names <- colnames(rna_data)
rna_data$time_group <- "未知"

# 使用细胞名称中的Day0和Day28前缀标识对照组和感染组
rna_data$time_group[grepl("^Day0_", cell_names)] <- "对照组"
rna_data$time_group[grepl("^Day28_", cell_names)] <- "感染组"

# 检查分组结果
print("时间点分组结果:")
print(table(rna_data$time_group))

# 如果分组结果不理想，尝试使用Groups变量
if(sum(rna_data$time_group %in% c("对照组", "感染组")) < 1000) {
  print("使用Groups变量进行分组...")
  # 查看Groups变量的值
  print("Groups变量的值:")
  print(table(rna_data$Groups))
  
  # 假设Groups变量包含day0和day28或类似值
  rna_data$time_group <- "未知"
  # 尝试常见的表示方法
  rna_data$time_group[grepl("day0|Day0|d0|D0|ctrl|control", rna_data$Groups, ignore.case = TRUE)] <- "对照组"
  rna_data$time_group[grepl("day28|Day28|d28|D28|infected|infection", rna_data$Groups, ignore.case = TRUE)] <- "感染组"
  
  print("使用Groups变量后的分组结果:")
  print(table(rna_data$time_group))
}

# 确保我们有足够的细胞被正确分组
if(sum(rna_data$time_group %in% c("对照组", "感染组")) < 1000) {
  stop("无法识别足够的对照组和感染组细胞，请手动检查数据")
}

# =============================================
# 绘制感染组vs对照组UMAP图 - 使用自定义颜色
# =============================================

# 提取UMAP坐标并添加组别信息
umap_data <- as.data.frame(rna_data[["umap"]]@cell.embeddings)
umap_data$group <- rna_data$time_group

# 只保留对照组和感染组的细胞
umap_data <- umap_data[umap_data$group %in% c("对照组", "感染组"), ]
print(paste0("用于绘图的细胞数: ", nrow(umap_data)))

# 设置自定义颜色 - 按照您的要求
group_colors <- c("对照组" = "#00bfc4", "感染组" = "#f7766c")

# 1. 基本UMAP图 - 按组别着色
p1 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = group)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  coord_fixed() +
  labs(
    title = "感染组 vs 对照组 细胞分布",
    color = "实验组别"
  )

# 保存图像
ggsave("感染vs对照组比较/感染vs对照组_UMAP.png", p1, width = 12, height = 10, dpi = 600)
ggsave("感染vs对照组比较/感染vs对照组_UMAP.pdf", p1, width = 12, height = 10)

# 2. 分开显示对照组和感染组
control_data <- umap_data[umap_data$group == "对照组", ]
infection_data <- umap_data[umap_data$group == "感染组", ]

print(paste0("对照组细胞数: ", nrow(control_data)))
print(paste0("感染组细胞数: ", nrow(infection_data)))

# 创建底图 - 所有细胞灰色背景
p_base <- ggplot() +
  geom_point(data = umap_data, aes(x = umap_1, y = umap_2), 
             color = "lightgrey", size = 0.5, alpha = 0.3) +
  theme_void() +
  coord_fixed() +
  theme(
    aspect.ratio = 1,
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# 对照组UMAP
p2 <- p_base +
  geom_point(data = control_data, aes(x = umap_1, y = umap_2), 
             color = group_colors["对照组"], size = 0.7, alpha = 0.8) +
  labs(title = "对照组细胞分布")

# 感染组UMAP
p3 <- p_base +
  geom_point(data = infection_data, aes(x = umap_1, y = umap_2), 
             color = group_colors["感染组"], size = 0.7, alpha = 0.8) +
  labs(title = "感染组细胞分布")

# 保存单独的图像
ggsave("感染vs对照组比较/对照组_UMAP.png", p2, width = 10, height = 8, dpi = 600)
ggsave("感染vs对照组比较/对照组_UMAP.pdf", p2, width = 10, height = 8)
ggsave("感染vs对照组比较/感染组_UMAP.png", p3, width = 10, height = 8, dpi = 600)
ggsave("感染vs对照组比较/感染组_UMAP.pdf", p3, width = 10, height = 8)

# 3. 组合图 - 并排比较
combined_plot <- p2 + p3 + plot_layout(ncol = 2)
ggsave("感染vs对照组比较/对照组vs感染组_并排UMAP.png", combined_plot, width = 16, height = 8, dpi = 600)
ggsave("感染vs对照组比较/对照组vs感染组_并排UMAP.pdf", combined_plot, width = 16, height = 8)

# 4. 两组细胞重叠显示，对比分布差异
p4 <- ggplot() +
  # 先绘制对照组
  geom_point(data = control_data, aes(x = umap_1, y = umap_2), 
             color = group_colors["对照组"], size = 0.7, alpha = 0.6) +
  # 再绘制感染组
  geom_point(data = infection_data, aes(x = umap_1, y = umap_2), 
             color = group_colors["感染组"], size = 0.7, alpha = 0.6) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  coord_fixed() +
  labs(title = "对照组与感染组细胞分布对比") +
  # 手动添加图例
  annotate("point", x = -Inf, y = -Inf, color = group_colors["对照组"], size = 3) +
  annotate("text", x = -Inf, y = -Inf, label = "对照组", hjust = -0.5, vjust = 0.5, size = 5) +
  annotate("point", x = -Inf, y = -Inf, color = group_colors["感染组"], size = 3) +
  annotate("text", x = -Inf, y = -Inf, label = "感染组", hjust = -0.5, vjust = 2.5, size = 5)

# 这个图未能正确添加图例，尝试另一种方法
# 创建带图例的版本
p4_legend <- ggplot() +
  geom_point(data = rbind(
    transform(control_data, type = "对照组"),
    transform(infection_data, type = "感染组")
  ), aes(x = umap_1, y = umap_2, color = type), size = 0.7, alpha = 0.6) +
  scale_color_manual(values = group_colors) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  coord_fixed() +
  labs(title = "对照组与感染组细胞分布对比", color = "实验组别")

ggsave("感染vs对照组比较/对照组vs感染组_重叠UMAP.png", p4_legend, width = 12, height = 10, dpi = 600)
ggsave("感染vs对照组比较/对照组vs感染组_重叠UMAP.pdf", p4_legend, width = 12, height = 10)

# 5. 使用细胞类型信息创建按组别分面的UMAP图
if("celltype_new" %in% colnames(rna_data@meta.data)) {
  # 添加细胞类型信息
  umap_data$celltype <- rna_data$celltype_new[match(rownames(umap_data), colnames(rna_data))]
  
  # 自定义颜色方案
  custom_colors <- c(
    "CD8 cells" = "#ef3f33",       # 红色
    "CD4 cells" = "#2f5885",       # 深蓝色
    "B cells" = "#26a870",         # 绿色
    "Plasma cells" = "#895da7",    # 紫色
    "NK cells" = "#3c9dd7",        # 天蓝色
    "DC cells" = "#00a18a",        # 青绿色
    "Macrophages" = "#f69886",     # 浅红色/粉橙色
    "Neutrophils" = "#d984b7",     # 粉紫色
    "Fibroblasts" = "#fdc182",     # 浅橙色
    "Endothelial cells" = "#a8d490", # 浅绿色
    "Epithelial cells" = "#d4c1ee"   # 浅紫色
  )
  
  # 按组别分面绘制UMAP图，着色为细胞类型
  p5 <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = celltype)) +
    geom_point(size = 0.7, alpha = 0.8) +
    scale_color_manual(values = custom_colors) +
    facet_wrap(~group, ncol = 2) +
    theme_void() +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = "感染组vs对照组细胞类型分布",
      color = "细胞类型"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  ggsave("感染vs对照组比较/感染vs对照组_细胞类型UMAP.png", p5, width = 16, height = 8, dpi = 600)
  ggsave("感染vs对照组比较/感染vs对照组_细胞类型UMAP.pdf", p5, width = 16, height = 8)
}

print("所有图像已保存在'感染vs对照组比较'文件夹中")