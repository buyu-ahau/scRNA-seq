# 设置工作目录
setwd("D:/onedrive/开题/结果/8、scRNA-seq/3.差异表达分析/基因符号_DE")

# 安装和加载必要的包
library(clusterProfiler)
library(org.Ss.eg.db)
library(ReactomePA)
library(ggplot2)
library(stringr)

# 创建结果目录
d <- "GO_KEGG_Results"
dir.create(d, showWarnings = FALSE)

# 读取所有细胞类型
cell_types <- c("B cells", "CD4 cells", "CD8 cells", "DC cells", "Endothelial cells",
                "Epithelial cells", "Fibroblasts", "Macrophages", "Neutrophils", "NK cells", "Plasma cells")

# 存储每个细胞类型的基因符号
gene_symbols_list <- list()

# 读取每个细胞类型的显著差异表达基因
for(cell in cell_types) {
  # 正确的文件路径
  file_path <- file.path(cell, paste0(cell, "_significant_DE_genes_with_symbols.csv"))
  
  if(file.exists(file_path)) {
    # 读取CSV文件
    genes_data <- read.csv(file_path)
    
    # 确保文件包含gene_symbol列
    if("gene_symbol" %in% colnames(genes_data)) {
      # 提取非空的基因符号（排除Ensembl ID格式的条目）
      gene_symbols <- genes_data$gene_symbol
      valid_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != "" & !grepl("^ENSSSCG", gene_symbols)]
      
      # 如果没有足够的有效基因符号，也尝试使用Ensembl ID
      if(length(valid_symbols) < 10) {
        cat(cell, "中没有足够的有效基因符号，尝试使用Ensembl ID\n")
        valid_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
      }
      
      # 移除了限制基因数量的代码
      
      gene_symbols_list[[cell]] <- valid_symbols
      cat("从", file_path, "中读取了", length(valid_symbols), "个有效基因符号\n")
    } else {
      warning(paste("文件", file_path, "中没有gene_symbol列。可用的列名有:", paste(colnames(genes_data), collapse=", ")))
    }
  } else {
    warning(paste("文件不存在:", file_path))
  }
}

# 打印每个细胞类型的有效基因数量
cat("\n每个细胞类型的基因数量:\n")
cell_gene_counts <- sapply(gene_symbols_list, length)
print(cell_gene_counts)

# 如果我们需要处理Ensembl ID，我们可以使用Ensembl_to_Symbol_mapping.csv
if(any(cell_gene_counts < 10) && file.exists("Ensembl_to_Symbol_mapping.csv")) {
  mapping_data <- read.csv("Ensembl_to_Symbol_mapping.csv")
  cat("\nEnsembl到Symbol映射文件的列名:", paste(colnames(mapping_data), collapse=", "), "\n")
  
  # 找出列名中包含"ensembl"和"symbol"的列
  ensembl_col <- NULL
  symbol_col <- NULL
  
  for(col in colnames(mapping_data)) {
    if(grepl("ensembl", col, ignore.case = TRUE)) ensembl_col <- col
    if(grepl("symbol", col, ignore.case = TRUE)) symbol_col <- col
  }
  
  if(!is.null(ensembl_col) && !is.null(symbol_col)) {
    # 创建映射字典
    ensembl_to_symbol <- setNames(mapping_data[[symbol_col]], mapping_data[[ensembl_col]])
    
    # 对有少量有效基因的细胞类型进行映射
    for(cell in names(gene_symbols_list)) {
      if(length(gene_symbols_list[[cell]]) < 10) {
        # 读取文件并获取所有Ensembl ID
        file_path <- file.path(cell, paste0(cell, "_significant_DE_genes_with_symbols.csv"))
        genes_data <- read.csv(file_path)
        ensembl_ids <- genes_data$gene
        
        # 映射Ensembl ID到基因符号
        symbols <- ensembl_to_symbol[ensembl_ids]
        valid_symbols <- symbols[!is.na(symbols) & symbols != ""]
        
        # 移除了限制基因数量的代码
        
        if(length(valid_symbols) >= 10) {
          gene_symbols_list[[cell]] <- valid_symbols
          cat("使用映射后，", cell, "现在有", length(valid_symbols), "个有效基因符号\n")
        }
      }
    }
  }
}

# 再次打印每个细胞类型的有效基因数量
cat("\n映射后每个细胞类型的基因数量:\n")
print(sapply(gene_symbols_list, length))

# 定义富集分析函数
com_go_kegg_ReactomePA_pig <- function(symbols_list, pro){ 
  # 转换基因符号为entrezID
  gcSample = lapply(symbols_list, function(y){ 
    if(length(y) == 0) return(character(0))
    
    # 尝试转换基因符号为EntrezID
    conversion <- tryCatch({
      AnnotationDbi::select(org.Ss.eg.db,
                            keys = y,
                            columns = "ENTREZID",
                            keytype = "SYMBOL")
    }, error = function(e) {
      message("转换基因符号时出错: ", e$message)
      return(data.frame(SYMBOL=character(0), ENTREZID=character(0)))
    })
    
    # 提取EntrezID，排除NA
    entrez_ids <- as.character(na.omit(conversion$ENTREZID))
    entrez_ids
  })
  
  # 检查转换后的基因数量
  cat("每个细胞类型转换后的基因数量:\n")
  entrez_counts <- sapply(gcSample, length)
  print(entrez_counts)
  
  # 移除空的细胞类型
  empty_cells <- entrez_counts == 0
  if(any(empty_cells)) {
    warning(paste("以下细胞类型没有有效的EntrezID，将被排除:", 
                  paste(names(gcSample)[empty_cells], collapse=", ")))
    gcSample <- gcSample[!empty_cells]
  }
  
  if(length(gcSample) == 0) {
    stop("没有有效的基因可以进行富集分析")
  }
  
  # KEGG通路分析
  tryCatch({
    xx <- compareCluster(gcSample, fun="enrichKEGG",
                         organism="ssc", pvalueCutoff=0.05)
    p1 <- dotplot(xx) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
      scale_y_discrete(labels=function(x) str_wrap(x, width=50))
    ggsave(file.path(d, paste0(pro,"_kegg.pdf")), plot=p1, width = 10, height = 8)
    write.csv(xx@compareClusterResult, file.path(d, paste0(pro,"_kegg.csv")))
    cat("KEGG分析完成\n")
  }, error = function(e) {
    message("KEGG分析出错: ", e$message)
  })
  
  # Reactome通路分析
  tryCatch({
    xx <- compareCluster(gcSample, fun="enrichPathway",
                         organism = "pig",
                         pvalueCutoff=0.05)
    p2 <- dotplot(xx) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
      scale_y_discrete(labels=function(x) str_wrap(x, width=50))
    ggsave(file.path(d, paste0(pro,"_ReactomePA.pdf")), plot=p2, width = 10, height = 8)
    write.csv(xx@compareClusterResult, file.path(d, paste0(pro,"_ReactomePA.csv")))
    cat("Reactome分析完成\n")
  }, error = function(e) {
    message("Reactome分析出错，可能不支持猪: ", e$message)
    # 尝试使用人类数据库
    tryCatch({
      message("尝试使用人类(human)进行Reactome分析...")
      xx <- compareCluster(gcSample, fun="enrichPathway",
                           organism = "human",
                           pvalueCutoff=0.05)
      p2 <- dotplot(xx) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
        scale_y_discrete(labels=function(x) str_wrap(x, width=50))
      ggsave(file.path(d, paste0(pro,"_ReactomePA_human.pdf")), plot=p2, width = 10, height = 8)
      write.csv(xx@compareClusterResult, file.path(d, paste0(pro,"_ReactomePA_human.csv")))
      cat("使用人类数据库的Reactome分析完成\n")
    }, error = function(e2) {
      message("使用人类数据库的Reactome分析也失败: ", e2$message)
    })
  })
  
  # GO生物过程(BP)分析
  tryCatch({
    formula_res <- compareCluster(
      gcSample, 
      fun="enrichGO", 
      OrgDb="org.Ss.eg.db",
      ont     = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
    
    lineage1_ego <- simplify(
      formula_res, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    p3 <- dotplot(lineage1_ego, showCategory=5) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
      scale_y_discrete(labels=function(x) str_wrap(x, width=50))
    ggsave(file.path(d, paste0(pro,"_GO_BP_cluster_simplified.pdf")), plot=p3, width = 15, height = 8)
    write.csv(lineage1_ego@compareClusterResult, file=file.path(d, paste0(pro,"_GO_BP_cluster_simplified.csv")))
    cat("GO BP分析完成\n")
  }, error = function(e) {
    message("GO BP分析出错: ", e$message)
  })
  
  # GO细胞组分(CC)分析
  tryCatch({
    formula_res <- compareCluster(
      gcSample, 
      fun="enrichGO", 
      OrgDb="org.Ss.eg.db",
      ont     = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
    
    lineage1_ego <- simplify(
      formula_res, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    p4 <- dotplot(lineage1_ego, showCategory=5) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
      scale_y_discrete(labels=function(x) str_wrap(x, width=50))
    ggsave(file.path(d, paste0(pro,"_GO_CC_cluster_simplified.pdf")), plot=p4, width = 15, height = 8)
    write.csv(lineage1_ego@compareClusterResult, file=file.path(d, paste0(pro,"_GO_CC_cluster_simplified.csv")))
    cat("GO CC分析完成\n")
  }, error = function(e) {
    message("GO CC分析出错: ", e$message)
  })
  
  # GO分子功能(MF)分析
  tryCatch({
    formula_res <- compareCluster(
      gcSample, 
      fun="enrichGO", 
      OrgDb="org.Ss.eg.db",
      ont     = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
    
    lineage1_ego <- simplify(
      formula_res, 
      cutoff=0.5, 
      by="p.adjust", 
      select_fun=min
    ) 
    p5 <- dotplot(lineage1_ego, showCategory=5) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
      scale_y_discrete(labels=function(x) str_wrap(x, width=50))
    ggsave(file.path(d, paste0(pro,"_GO_MF_cluster_simplified.pdf")), plot=p5, width = 15, height = 8)
    write.csv(lineage1_ego@compareClusterResult, file=file.path(d, paste0(pro,"_GO_MF_cluster_simplified.csv")))
    cat("GO MF分析完成\n")
  }, error = function(e) {
    message("GO MF分析出错: ", e$message)
  })
  
  cat("富集分析完成。结果保存在", d, "目录中。\n")
}

# 运行富集分析
com_go_kegg_ReactomePA_pig(gene_symbols_list, pro="Pig_CellTypes")
