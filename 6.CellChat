##################################################################
# CellChat分析完整流程代码
# 作者: plmmmz
# 日期: 2025-08-20
##################################################################

# 1. 基本设置与加载包
#-----------------------------------------------------------------
library(CellChat)
library(patchwork)
library(ggplot2)
library(qs)
library(Seurat)

# 设置工作目录（如需要）
# setwd("/disk192/users_dir/buyu/1.布宇/3.布宇scATAC-seq/2025.8.14/3.CellChat/")

# 2. 数据加载与预处理
#-----------------------------------------------------------------
# 加载转换后的人类同源基因数据
rna_human <- qs::qread("scRNA_pig_human.qs")
print("数据加载完成")

# 准备数据
data.input <- GetAssayData(rna_human, slot = "data", assay = "RNA")  # 使用规范化数据
meta.data <- rna_human@meta.data

# 创建CellChat对象
cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "celltype")
print("CellChat对象创建完成")

# 查看细胞类型分布
groupSize <- as.numeric(table(cellchat@idents))
print("细胞类型分布:")
print(table(cellchat@idents))

# 3. 加载数据库并预处理
#-----------------------------------------------------------------
# 加载人类CellChatDB数据库
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# 设置使用的数据库
cellchat@DB <- CellChatDB

# 预处理数据
cellchat <- subsetData(cellchat)
print(paste("从", nrow(data.input), "个基因中筛选出", 
            nrow(cellchat@data.signaling), "个在CellChatDB中的信号基因"))

# 识别过表达的配体-受体对
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
print(paste("用于信号推断的高变异配体-受体对数量:", 
            nrow(cellchat@LR$LRsig)))

# 4. 计算细胞通讯
#-----------------------------------------------------------------
# 计算通讯概率
cellchat <- computeCommunProb(cellchat)
print("通讯概率计算完成")

# 过滤低质量的细胞-细胞通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 计算通路级别的通讯概率
cellchat <- computeCommunProbPathway(cellchat)
print("通路级别细胞通讯推断完成")

# 聚合细胞通讯网络
cellchat <- aggregateNet(cellchat)
print("聚合细胞通讯网络计算完成")

# 保存CellChat对象
qs::qsave(cellchat, "cellchat_pig_human_results.qs")
print("CellChat分析结果已保存到: cellchat_pig_human_results.qs")

# 5. 可视化分析: 整体网络
#-----------------------------------------------------------------
# 创建保存图像的目录
dir.create("cellchat_plots", showWarnings = FALSE)

# 可视化总体通讯网络
pdf("cellchat_plots/1_overall_network.pdf", width = 12, height = 10)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()

# 可视化每种细胞类型的交互数量
pdf("cellchat_plots/2_interaction_counts.pdf", width = 10, height = 8)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                    edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# 6. 可视化分析: 信号模式
#-----------------------------------------------------------------
# 计算网络中心度
cellchat <- netAnalysis_computeCentrality(cellchat)
print("网络中心度计算完成")

# 热图显示细胞的信号发送/接收模式
pdf("cellchat_plots/3_signaling_role.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(cellchat, width = 8, height = 6, font.size = 10)
dev.off()

# 热图展示细胞发送和接收的模式
pdf("cellchat_plots/4_signaling_role_heatmap.pdf", width = 15, height = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 10, title = "Outgoing signaling")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 10, title = "Incoming signaling")
ht1 + ht2
dev.off()

# 7. 可视化分析: 通路重要性
#-----------------------------------------------------------------
# 获取信号通路并排序
pathways <- cellchat@netP$pathways
print("可用的信号通路:")
print(pathways)

# 按照通路重要性排序
pathway_importance <- sapply(pathways, function(x) {
  sum(cellchat@netP$weight[[x]])
})
pathways_ordered <- names(sort(pathway_importance, decreasing = TRUE))
top_pathways <- pathways_ordered[1:min(10, length(pathways_ordered))]
print("前10个重要通路:")
print(top_pathways)

# 可视化前10个重要通路 - 替代方法（网络图）
pdf("cellchat_plots/5_top_pathways_alt.pdf", width = 12, height = 8)
for (i in 1:length(top_pathways)) {
  pathway <- top_pathways[i]
  netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  title(main = pathway)
}
dev.off()

# 可视化每个通路的网络
pdf("cellchat_plots/6_pathway_networks.pdf", width = 12, height = 9)
par(mfrow = c(1,1))
for (pathway in top_pathways) {
  netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  title(main = pathway, line = -1.5)
}
dev.off()

# 可视化每个通路的发送-接收关系
pdf("cellchat_plots/7_pathway_contribution.pdf", width = 12, height = 8)
for (pathway in top_pathways) {
  print(paste("生成通路热图:", pathway))
  netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
}
dev.off()

# 8. 可视化分析: 配体-受体对
#-----------------------------------------------------------------
# 提取所有配体-受体对信息
lr_pairs <- subsetCommunication(cellchat)
# 按照权重排序
lr_pairs <- lr_pairs[order(lr_pairs$prob, decreasing = TRUE),]
print("前20个重要的配体-受体对:")
print(head(lr_pairs, 20))

# 可视化前5个配体-受体对的通讯模式
pdf("cellchat_plots/8_top_lr_pairs.pdf", width = 12, height = 10)
top_pairs <- head(lr_pairs, 5)
for (i in 1:nrow(top_pairs)) {
  pair <- top_pairs[i, ]
  netVisual_individual(cellchat, signaling = pair$pathway_name, 
                      pairLR.use = paste0(pair$ligand, "_", pair$receptor))
}
dev.off()

# 9. 可视化分析: 细胞通讯模式
#-----------------------------------------------------------------
# 创建通路交流模式的热图（无需UMAP）
pdf("cellchat_plots/9_pathway_patterns.pdf", width = 12, height = 10)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", title = "Overall signaling patterns")
dev.off()

# 创建气泡图显示通路重要性
pdf("cellchat_plots/10_pathway_importance_bubble.pdf", width = 12, height = 9)
for (i in 1:min(5, length(top_pathways))) {
  pathway <- top_pathways[i]
  netVisual_bubble(cellchat, sources.use = 1:nrow(cellchat@net$weight), 
                  targets.use = 1:ncol(cellchat@net$weight), signaling = pathway, 
                  title = paste0("Signaling pathway: ", pathway))
}
dev.off()

# 10. 可视化分析: 特定通路的细胞互作
#-----------------------------------------------------------------
# 和弦图显示特定通路下的细胞互作
pdf("cellchat_plots/11_detailed_interactions.pdf", width = 12, height = 10)
for (i in 1:min(3, length(top_pathways))) {
  pathway <- top_pathways[i]
  netVisual_chord_cell(cellchat, signaling = pathway, title = paste0(pathway, " signaling network"))
}
dev.off()

# 获取细胞类型列表
cell_types <- levels(cellchat@idents)
print("细胞类型:")
print(cell_types)

# 免疫细胞之间的互作
immune_cells <- c("Macrophages", "CD4 cells", "CD8 cells", "B cells", "NK cells")
immune_idx <- which(cell_types %in% immune_cells)

pdf("cellchat_plots/12_immune_cell_interactions.pdf", width = 12, height = 10)
if (length(immune_idx) >= 2) {
  for (i in 1:min(3, length(top_pathways))) {
    pathway <- top_pathways[i]
    tryCatch({
      netVisual_chord_cell(cellchat, signaling = pathway, 
                         sources.use = immune_idx, targets.use = immune_idx,
                         title.name = paste0(pathway, " signaling in immune cells"))
    }, error = function(e) {
      print(paste("无法为", pathway, "生成图表:", e$message))
    })
  }
}
dev.off()

# 结构细胞之间的互作
structural_cells <- c("Fibroblasts", "Endothelial cells", "Epithelial cells")
structural_idx <- which(cell_types %in% structural_cells)

pdf("cellchat_plots/13_structural_cell_interactions.pdf", width = 12, height = 10)
if (length(structural_idx) >= 2) {
  for (i in 1:min(3, length(top_pathways))) {
    pathway <- top_pathways[i]
    tryCatch({
      netVisual_chord_cell(cellchat, signaling = pathway, 
                         sources.use = structural_idx, targets.use = structural_idx,
                         title.name = paste0(pathway, " signaling in structural cells"))
    }, error = function(e) {
      print(paste("无法为", pathway, "生成图表:", e$message))
    })
  }
}
dev.off()

# 11. 结果导出
#-----------------------------------------------------------------
# 创建结果文件夹
dir.create("cellchat_results", showWarnings = FALSE)

# 导出重要的配体-受体对
write.csv(lr_pairs, "cellchat_results/important_LR_pairs.csv", row.names = FALSE)

# 导出细胞类型之间的通讯强度矩阵
write.csv(cellchat@net$weight, "cellchat_results/cell_communication_strength.csv")

# 导出通路的重要性得分
pathway_scores <- data.frame(
  Pathway = names(pathway_importance),
  Importance = pathway_importance,
  stringsAsFactors = FALSE
)
pathway_scores <- pathway_scores[order(pathway_scores$Importance, decreasing = TRUE),]
write.csv(pathway_scores, "cellchat_results/pathway_importance.csv", row.names = FALSE)

# 导出细胞类型的信号角色信息
cell_roles <- data.frame(
  Cell_Type = rownames(cellchat@netP$centr),
  Outgoing = cellchat@netP$centr$outdeg,
  Incoming = cellchat@netP$centr$indeg,
  Influence = cellchat@netP$centr$influence,
  stringsAsFactors = FALSE
)
write.csv(cell_roles, "cellchat_results/cell_signaling_roles.csv", row.names = FALSE)

# 12. 汇总可视化
#-----------------------------------------------------------------
# 散点图展示发送和接收模式
pdf("cellchat_plots/16_signaling_roles_scatter.pdf", width = 10, height = 8)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# 总体通讯强度热图
pdf("cellchat_plots/17_overall_communication_heatmap.pdf", width = 12, height = 10)
netVisual_heatmap(cellchat, measure = "weight", title = "Overall communication strength")
dev.off()

# 热图显示每个通路的贡献
pdf("cellchat_plots/18_pathway_contribution_heatmap.pdf", width = 15, height = 12)
ht <- netAnalysis_contribution(cellchat, signaling = top_pathways)
dev.off()

# 13. 保存完整分析结果
#-----------------------------------------------------------------
qs::qsave(cellchat, "cellchat_pig_human_complete.qs")

print("分析完成！结果已保存到cellchat_plots和cellchat_results文件夹")
