# 1. 加载必要的库
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pathview)
library(dplyr)
# 设置工作目录（根据你的路径修改）
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")
# ---------------------------
# 步骤 2：数据读取与筛选
# ---------------------------

# 定义一个自动化筛选函数
get_sig_genes <- function(file_path, platform = "TMT") {
  df <- read.csv(file_path)
  
  # 根据平台设定初步的 logFC 阈值 (TMT压缩效应更明显)
  lfc_threshold <- if(platform == "TMT") 0.25 else 0.58 
  
  # 筛选 adj.P.Val < 0.05 的差异蛋白
  sig_df <- df %>% 
    filter(adj.P.Val < 0.05 & abs(logFC) > lfc_threshold)
  
  cat(paste0(platform, " 筛选出的差异蛋白数量: ", nrow(sig_df), "\n"))
  return(sig_df)
}

# 运行筛选 
tmt_sig <- get_sig_genes("tmt_EMT_vs_Epithelial_Corrected.csv", "TMT")
lf_sig <- get_sig_genes("LF_EMT_vs_Epithelial_Corrected.csv", "LF")

# ---------------------------
# 步骤 3：ID 转换 (Symbol -> Entrez ID)
# ---------------------------

convert_id <- function(sig_df) {
  gene_map <- bitr(sig_df$Protein, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  merged <- inner_join(sig_df, gene_map, by = c("Protein" = "SYMBOL"))
  return(merged)
}

tmt_mapped <- convert_id(tmt_sig)
lf_mapped <- convert_id(lf_sig)

# ---------------------------
# 步骤 4：KEGG 富集分析与气泡图
# ---------------------------

run_kegg <- function(mapped_df, title_name) {
  kk <- enrichKEGG(gene = mapped_df$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
  
  # 绘制气泡图
  p <- dotplot(kk, showCategory = 20) + 
    ggtitle(paste("KEGG Enrichment -", title_name)) +
    theme_bw()
  
  print(p)
  return(kk)
}
# 在运行任何 enrichKEGG 之前运行这一行
options(timeout = 300)
kk_tmt <- run_kegg(tmt_mapped, "TMT Platform")
kk_lf <- run_kegg(lf_mapped, "LF Platform")

# ---------------------------
# 步骤 5：Pathview 通路可视化
# ---------------------------

# 准备包含 logFC 的命名向量
genelist_tmt <- tmt_mapped$logFC
names(genelist_tmt) <- tmt_mapped$ENTREZID

# 选取 TMT 结果中最显著的第一个通路 ID 进行绘制
# 例如：hsa04510 (Focal adhesion)
top_pathway <- kk_tmt@result$ID[1] 

pathview(gene.data = genelist_tmt, 
         pathway.id = top_pathway, 
         species = "hsa",
         out.suffix = "TMT_Pathview",
         low = list(gene = "green"), 
         high = list(gene = "red"))

# ---------------------------
# 步骤 2：数据读取与筛选
# 运行筛选 Hyper_vs_EMT
tmt_sig1 <- get_sig_genes("tmt_Hyper_vs_EMT_Corrected.csv", "TMT")
lf_sig1 <- get_sig_genes("LF_Hyper_vs_EMT_Corrected.csv", "LF")
# ---------------------------
# 步骤 3：ID 转换 (Symbol -> Entrez ID)
# ---------------------------

convert_id <- function(sig_df) {
  gene_map <- bitr(sig_df$Protein, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  merged <- inner_join(sig_df, gene_map, by = c("Protein" = "SYMBOL"))
  return(merged)
}

tmt_mapped1 <- convert_id(tmt_sig1)
lf_mapped1 <- convert_id(lf_sig1)

# ---------------------------
# 步骤 4：KEGG 富集分析与气泡图
# ---------------------------

run_kegg <- function(mapped_df, title_name) {
  kk <- enrichKEGG(gene = mapped_df$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
  
  # 绘制气泡图
  p <- dotplot(kk, showCategory = 20) + 
    ggtitle(paste("KEGG Enrichment -", title_name)) +
    theme_bw()
  
  print(p)
  return(kk)
}
# 在运行任何 enrichKEGG 之前运行这一行
options(timeout = 300)
kk_tmt1 <- run_kegg(tmt_mapped1, "TMT Platform")
kk_lf1 <- run_kegg(lf_mapped1, "LF Platform")

# ---------------------------
# 步骤 2：数据读取与筛选
# 运行筛选 Hyper_vs_EPi
tmt_sig2 <- get_sig_genes("tmt_Hyper_vs_Epithelial_Corrected.csv", "TMT")
lf_sig2 <- get_sig_genes("LF_Hyper_vs_Epithelial_Corrected.csv", "LF")
# ---------------------------
# 步骤 3：ID 转换 (Symbol -> Entrez ID)
# ---------------------------

convert_id <- function(sig_df) {
  gene_map <- bitr(sig_df$Protein, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  merged <- inner_join(sig_df, gene_map, by = c("Protein" = "SYMBOL"))
  return(merged)
}

tmt_mapped2 <- convert_id(tmt_sig2)
lf_mapped2 <- convert_id(lf_sig2)

# ---------------------------
# 步骤 4：KEGG 富集分析与气泡图
# ---------------------------

run_kegg <- function(mapped_df, title_name) {
  kk <- enrichKEGG(gene = mapped_df$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
  
  # 绘制气泡图
  p <- dotplot(kk, showCategory = 20) + 
    ggtitle(paste("KEGG Enrichment -", title_name)) +
    theme_bw()
  
  print(p)
  return(kk)
}
# 在运行任何 enrichKEGG 之前运行这一行
options(timeout = 300)
kk_tmt2 <- run_kegg(tmt_mapped2, "TMT Platform")
kk_lf2 <- run_kegg(lf_mapped2, "LF Platform")

# 加载必要的包
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db) # 假设物种为人类
library(dplyr)

# 读取 TMT 和 LF 的差异分析结果 EMT vs Epi
tmt_data <- read.csv("tmt_EMT_vs_Epithelial_Corrected.csv")
lf_data <- read.csv("LF_EMT_vs_Epithelial_Corrected.csv")

# --- 步骤 1: ID 转换 (Symbol -> Entrez ID) ---
# TMT 数据转换
tmt_ids <- bitr(tmt_data$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
tmt_merged <- merge(tmt_data, tmt_ids, by.x="Protein", by.y="SYMBOL")

# LF 数据转换
lf_ids <- bitr(lf_data$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
lf_merged <- merge(lf_data, lf_ids, by.x="Protein", by.y="SYMBOL")

# --- 步骤 2: 合并两个平台的数据 ---
# 仅保留 EntrezID 和 logFC 列
tmt_sub <- tmt_merged[, c("ENTREZID", "logFC")]
colnames(tmt_sub)[2] <- "TMT_logFC"

lf_sub <- lf_merged[, c("ENTREZID", "logFC")]
colnames(lf_sub)[2] <- "LF_logFC"

# 使用 full_join 合并，确保两个平台检测到的蛋白都能涵盖
combined_df <- full_join(tmt_sub, lf_sub, by="ENTREZID")

# --- 步骤 3: 构建 Pathview 矩阵 ---
# 转换为矩阵，行名为 EntrezID
pv_matrix <- as.matrix(combined_df[, c("TMT_logFC", "LF_logFC")])
rownames(pv_matrix) <- combined_df$ENTREZID

# 处理 NA 值（可选：若某一平台未检测到，设为 0 即灰色，或保持 NA）
pv_matrix[is.na(pv_matrix)] <- 0
pathview(gene.data  = pv_matrix, 
         pathway.id = "hsa04810", 
         species    = "hsa",
         kegg.native = TRUE,     # 生成经典的 KEGG 风格图
         same.layer  = TRUE,     # 确保颜色条和图层一致
         limit      = list(gene = 2, cpd = 1), # 设置颜色映射范围 (logFC: -2 到 2)
         low        = "green",   # 下调颜色
         mid        = "gray",    # 无差异颜色
         high       = "red",     # 上调颜色
         out.suffix = "TMT_vs_LF", # 输出文件名后缀
         multi.state = TRUE      # 关键参数：支持多列对比
)
# 1. 加载库
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# 2. 读取你的数据
tmt_df1 <- read.csv("tmt_EMT_vs_Epithelial_Corrected.csv")
lf_df1 <- read.csv("LF_EMT_vs_Epithelial_Corrected.csv")

# 3. ID 转换 (Symbol -> Entrez ID)
tmt_ids1 <- bitr(tmt_df1$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
lf_ids1 <- bitr(lf_df1$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# 4. 合并数据并提取 logFC
# TMT 部分
tmt_data1 <- merge(tmt_df1, tmt_ids1, by.x="Protein", by.y="SYMBOL") %>%
  dplyr::select(ENTREZID, logFC_TMT = logFC)

# LF 部分
lf_data1 <- merge(lf_df1, lf_ids1, by.x="Protein", by.y="SYMBOL") %>%
  dplyr::select(ENTREZID, logFC_LF = logFC)

# 全连接合并
combined <- full_join(tmt_data1, lf_data1, by="ENTREZID")

# 5. 构建 Pathview 矩阵 (行名为 EntrezID，列为不同平台)
pv_matrix1 <- as.matrix(combined[, 2:3])
rownames(pv_matrix1) <- combined$ENTREZID

# 6. 生成可视化图 (Regulation of actin cytoskeleton)
pathview(gene.data  = pv_matrix1, 
         pathway.id = "hsa04810", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF EMT VS EPI 1")
# 6. 生成可视化图 (Cytoskeleton in muscle cells)
pathview(gene.data  = pv_matrix1, 
         pathway.id = "hsa05410", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF EMT VS EPI 2")

#调用 Pathview 绘图（ECM）
pathview(gene.data  = pv_matrix1, 
         pathway.id = "hsa04512", # ECM-receptor interaction
         species    = "hsa",
         multi.state = TRUE,      # 启用多列对比显示
         same.layer  = TRUE,      # 保持图层一致
         limit      = list(gene = 2, cpd = 1), # 设置 logFC 范围为 -2 到 2
         bins       = list(gene = 20),         # 增加颜色梯度细腻度
         low        = "green",    # 下调
         mid        = "gray",     # 无显著差异
         high       = "red",      # 上调
         out.suffix = "TMT_vs_LF_EMT VS EPI 3"
)

# 1. 加载库
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
options(timeout = 600)
# 2. 读取你的数据
tmt_df2 <- read.csv("tmt_Hyper_vs_EMT_Corrected.csv")
lf_df2 <- read.csv("LF_Hyper_vs_EMT_Corrected.csv")

# 3. ID 转换 (Symbol -> Entrez ID)
tmt_ids2 <- bitr(tmt_df2$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
lf_ids2 <- bitr(lf_df2$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# 4. 合并数据并提取 logFC
# TMT 部分
tmt_data2 <- merge(tmt_df2, tmt_ids2, by.x="Protein", by.y="SYMBOL") %>%
  dplyr::select(ENTREZID, logFC_TMT = logFC)

# LF 部分
lf_data2 <- merge(lf_df2, lf_ids2, by.x="Protein", by.y="SYMBOL") %>%
  dplyr::select(ENTREZID, logFC_LF = logFC)

# 全连接合并
combined2 <- full_join(tmt_data2, lf_data2, by="ENTREZID")

# 5. 构建 Pathview 矩阵 (行名为 EntrezID，列为不同平台)
pv_matrix2 <- as.matrix(combined2[, 2:3])
rownames(pv_matrix2) <- combined2$ENTREZID

# 6. 生成可视化图 (Hypertrophic cardiomyopathy)
pathview(gene.data  = pv_matrix2, 
         pathway.id = "hsa05414", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF_hyper vs EMT2")
# 6. 生成可视化图 (DNA replication)
pathview(gene.data  = pv_matrix2, 
         pathway.id = "hsa03030", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF_hyper vs EMT3")

# 2. 读取你的数据hsa00071
tmt_df3 <- read.csv("tmt_Hyper_vs_Epithelial_Corrected.csv")
lf_df3 <- read.csv("LF_Hyper_vs_Epithelial_Corrected.csv")

# 3. ID 转换 (Symbol -> Entrez ID)
tmt_ids3 <- bitr(tmt_df3$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
lf_ids3 <- bitr(lf_df3$Protein, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# 4. 合并数据并提取 logFC
# TMT 部分
tmt_data3 <- merge(tmt_df3, tmt_ids3, by.x="Protein", by.y="SYMBOL") %>%
  dplyr::select(ENTREZID, logFC_TMT = logFC)

# LF 部分
lf_data3 <- merge(lf_df3, lf_ids3, by.x="Protein", by.y="SYMBOL") %>%
  dplyr::select(ENTREZID, logFC_LF = logFC)

# 全连接合并
combined3 <- full_join(tmt_data3, lf_data3, by="ENTREZID")

# 5. 构建 Pathview 矩阵 (行名为 EntrezID，列为不同平台)
pv_matrix3 <- as.matrix(combined3[, 2:3])
rownames(pv_matrix3) <- combined3$ENTREZID
colnames(combined3)
# 确保只包含数字，且行名是字符型的 EntrezID
pv_matrix3 <- as.matrix(combined3[, c("logFC_TMT", "logFC_LF")])
rownames(pv_matrix3) <- as.character(combined3$ENTREZID)
# 关键检查：确保没有任何 NA
if(any(is.na(pv_matrix3))) stop("矩阵中仍存在 NA 值！")
# 6. 生成可视化图 (Fatty acid metabolism)
pathview(gene.data  = pv_matrix3, 
         pathway.id = "hsa00071", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF_hyper vs Epi 1")

# 6. 生成可视化图 (carbon metabolism)
pathview(gene.data  = pv_matrix3, 
         pathway.id = "hsa01200", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         match.data = FALSE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF_hyper vs Epi 2")

# 6. 生成可视化图 (carbon metabolism)
pathview(gene.data  = pv_matrix2, 
         pathway.id = "hsa01100", 
         species    = "hsa",
         multi.state = TRUE,    # 开启多状态对比模式
         same.layer  = TRUE,
         limit       = list(gene = 1.5, cpd = 1), # 颜色范围 logFC (-1.5 to 1.5)
         low         = "green", 
         mid         = "gray", 
         high        = "red",
         out.suffix  = "TMT_vs_LF_hyper vs Epi 3")
