# 猪scRNA-seq数据CellChat分析报告

**作者:** plmmmz  
**日期:** 2025-08-20

## 分析概述

本报告展示了对猪scRNA-seq数据（转换为人类同源基因）进行的CellChat细胞通讯分析结果。通过这一分析，我们鉴定了不同细胞类型之间的交互模式、重要的信号通路以及关键的配体-受体对。

### 数据统计
- **细胞类型数量:** 11
- **识别的通路数量:** 37
- **配体-受体对数量:** 476
- **生成的图表数量:** 20

## 主要细胞类型

分析中包含的细胞类型如下：

| 细胞类型 | 数量 |
|---------|------|
| Macrophages | - |
| CD4 cells | - |
| CD8 cells | - |
| B cells | - |
| Epithelial cells | - |
| Fibroblasts | - |
| Endothelial cells | - |
| Neutrophils | - |
| Plasma cells | - |
| NK cells | - |
| DC cells | - |

## 关键信号通路

通过计算通路的累积通讯强度，我们识别出以下前10个重要通路：

| 排名 | 通路 | 重要性得分 |
|-----|------|-----------|
| 1 | COLLAGEN | - |
| 2 | LAMININ | - |
| 3 | MK | - |
| 4 | FN1 | - |
| 5 | APP | - |
| 6 | TGFb | - |
| 7 | ADGRE | - |
| 8 | PTPRM | - |
| 9 | MHC-II | - |
| 10 | PECAM1 | - |

## 重要的配体-受体对

根据通讯强度排序的前10个重要配体-受体对：

| 来源细胞 | 目标细胞 | 配体 | 受体 | 通路 | 通讯强度 |
|---------|---------|-----|------|------|---------|
| Fibroblasts | Macrophages | COL1A1 | CD44 | COLLAGEN | 0.355 |
| Endothelial cells | Endothelial cells | PECAM1 | PECAM1 | PECAM1 | 0.335 |
| Endothelial cells | Endothelial cells | PTPRM | PTPRM | PTPRM | 0.331 |
| Fibroblasts | CD4 cells | COL1A1 | CD44 | COLLAGEN | 0.284 |
| Fibroblasts | Fibroblasts | COL1A1 | ITGA9_ITGB1 | COLLAGEN | 0.271 |
| Fibroblasts | NK cells | COL1A1 | CD44 | COLLAGEN | 0.268 |
| Fibroblasts | CD8 cells | COL1A1 | CD44 | COLLAGEN | 0.266 |
| Fibroblasts | Plasma cells | COL1A1 | CD44 | COLLAGEN | 0.265 |
| Fibroblasts | DC cells | COL1A1 | CD44 | COLLAGEN | 0.265 |
| Fibroblasts | Macrophages | COL4A2 | CD44 | COLLAGEN | 0.263 |

## 细胞通讯角色

细胞类型作为信号发送者和接收者的角色：

| 细胞类型 | 传出信号 | 传入信号 | 影响力 |
|---------|---------|---------|-------|
| Fibroblasts | - | - | - |
| Endothelial cells | - | - | - |
| Macrophages | - | - | - |
| CD4 cells | - | - | - |
| CD8 cells | - | - | - |
| B cells | - | - | - |
| Epithelial cells | - | - | - |
| Neutrophils | - | - | - |
| NK cells | - | - | - |
| Plasma cells | - | - | - |
| DC cells | - | - | - |

## 主要研究发现

### 关键通路

1. **胶原蛋白(COLLAGEN)通路**: 最重要的通讯通路，主要由成纤维细胞发出。这表明细胞外基质组分在组织中的细胞通讯中发挥着核心作用。
2. **层粘连蛋白(LAMININ)通路**: 第二重要的通路，同样与细胞外基质相关，参与细胞黏附和组织架构的维持。
3. **中间激酶(MK)通路**: 在细胞增殖、分化和存活中具有重要作用。
4. **纤连蛋白(FN1)通路**: 作为细胞外基质蛋白，在细胞黏附、生长和迁移中发挥重要作用。
5. **APP通路**: 淀粉样前体蛋白相关信号，可能与组织稳态和细胞功能调节有关。
6. **TGFb通路**: 转化生长因子β信号通路，参与细胞增殖、分化和免疫调节。

### 重要的配体-受体对

1. **COL1A1-CD44**: 成纤维细胞向巨噬细胞发送的最强信号，表明成纤维细胞产生的I型胶原与多种免疫细胞上的CD44受体有重要互作。
2. **PECAM1-PECAM1**: 内皮细胞之间的自分泌信号，与血管稳态和内皮完整性维持相关。
3. **PTPRM-PTPRM**: 内皮细胞之间的另一重要自分泌信号，参与细胞-细胞接触和细胞黏附。
4. **COL1A1-ITGA9_ITGB1**: 成纤维细胞自分泌通讯，可能参与细胞外基质重塑和组织修复。
5. **APP-SORL1**: 内皮细胞与中性粒细胞之间的通讯，在免疫反应和组织稳态中可能发挥作用。

### 关键细胞类型

1. **成纤维细胞(Fibroblasts)**: 主要的信号发送者，特别是通过胶原蛋白通路，向多种细胞类型发送信号，在组织微环境中起核心调控作用。
2. **巨噬细胞(Macrophages)**: 重要的信号接收者，接收来自成纤维细胞的多种胶原蛋白信号，提示间质-免疫细胞互作在组织中的重要性。
3. **内皮细胞(Endothelial cells)**: 既发送又接收信号，具有强烈的自分泌通讯（PECAM1和PTPRM），表明血管内皮细胞具有自我调节能力。
4. **CD4和CD8 T细胞**: 作为重要的信号接收者，尤其接收来自成纤维细胞的胶原蛋白信号。

## 可视化结果

分析产生了20个可视化PDF文件，存放在`cellchat_plots`文件夹中，主要包括：

### 整体网络可视化
- **1_overall_network.pdf**: 总体细胞通讯网络圆形布局图，显示互作数量和强度
- **2_interaction_counts.pdf**: 每种细胞类型的发送信号分布图
- **17_overall_communication_heatmap.pdf**: 总体通讯强度热图

### 细胞角色分析
- **3_signaling_role.pdf**: 细胞类型的信号角色网络图
- **4_signaling_role_heatmap.pdf**: 细胞发送和接收信号的热图
- **9_pathway_patterns.pdf**: 所有通路的信号模式热图
- **16_signaling_roles_scatter.pdf**: 发送者和接收者模式散点图

### 通路重要性分析
- **5_top_pathways_alt.pdf**: 前10个重要通路的网络图
- **6_pathway_networks.pdf**: 每个重要通路的详细网络图
- **7_pathway_contribution.pdf**: 每个通路的发送-接收关系热图
- **10_pathway_importance_bubble.pdf**: 通路重要性气泡图
- **18_pathway_contribution_heatmap.pdf**: 通路贡献热图

### 特定细胞类型分析
- **11_detailed_interactions.pdf**: 特定通路的细胞互作和弦图
- **12_immune_cell_interactions.pdf**: 免疫细胞之间的互作分析
- **13_structural_cell_interactions.pdf**: 结构细胞之间的互作分析

## 生物学解释与意义

1. **细胞外基质信号的主导地位**: 
   - COLLAGEN、LAMININ和FN1通路的重要性表明细胞外基质在组织中起关键调控作用
   - 成纤维细胞分泌的胶原蛋白通过CD44受体调控多种免疫细胞功能，形成"基质-免疫轴"
   - 这种通讯可能参与组织修复、炎症调控和免疫反应

2. **免疫-间质细胞相互作用**:
   - 成纤维细胞和巨噬细胞之间的强互作提示组织修复或炎症调控的重要机制
   - CD44作为多种胶原蛋白的受体，在免疫细胞中扮演关键角色
   - T细胞（CD4和CD8）也受到成纤维细胞信号的强烈影响，提示适应性免疫反应可能受基质细胞调控

3. **内皮细胞的自分泌调节**:
   - PECAM1和PTPRM的自分泌信号表明内皮细胞具有自我调节能力
   - 这种自分泌信号对于维持血管完整性、内皮细胞连接和血管稳态可能至关重要
   - 同时，内皮细胞也通过APP-SORL1等信号与中性粒细胞等免疫细胞通讯

4. **TGFb信号通路的重要性**:
   - 作为第六大重要通路，TGFb参与多种生物学过程，包括细胞增殖、分化和免疫调节
   - 在组织修复、纤维化和免疫抑制中可能发挥关键作用

## 结论与展望

本分析揭示了猪组织中复杂的细胞通讯网络，特别是成纤维细胞、巨噬细胞和内皮细胞之间的重要互作。细胞外基质相关通路（如胶原蛋白和层粘连蛋白）在细胞通讯中起主导作用，提示它们在组织稳态和功能中的关键地位。

未来研究可以：
1. 通过体外实验验证关键的配体-受体对（如COL1A1-CD44）在特定细胞类型中的功能
2. 探索特定疾病状态下这些细胞通讯网络的变化
3. 利用药物干预特定通路（如TGFb）来调控细胞通讯，可能为疾病治疗提供新策略
4. 比较猪与人类细胞通讯网络的异同，验证猪模型在人类疾病研究中的应用价值

这些发现为理解组织微环境中的细胞间相互作用提供了宝贵见解，有助于后续功能验证和潜在干预策略的设计。

## 数据和代码可用性

- 完整的分析结果保存在`cellchat_results`文件夹中
- 所有可视化图表保存在`cellchat_plots`文件夹中
- 完整的CellChat对象保存为`cellchat_pig_human_complete.qs`，可用于后续分析
