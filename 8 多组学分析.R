#############@转录组和蛋白组共分析#####################
library(limma)
library(tidyverse)
# 设置工作目录（根据你的路径修改）
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")
# 通用差异分析函数
get_degs <- function(file_path, meta_data, label) {
  # 读取数据 (自动识别 gz 压缩)
  df <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  # 对齐样本
  common <- intersect(colnames(df), meta_data$sample_id)
  df_sub <- df[, common]
  m_sub <- meta_data[match(common, meta_data$sample_id), ]
  
  # Limma 差异分析
  grp <- factor(m_sub$phenotype, levels = c("Epithelial", "EMT", "Hypermutated"))
  des <- model.matrix(~0 + grp)
  colnames(des) <- levels(grp)
  
  fit <- lmFit(df_sub, des)
  cont <- makeContrasts(EMT_vs_Epi = EMT - Epithelial, levels = des)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  
  # 返回结果
  res <- topTable(fit2, coef = "EMT_vs_Epi", number = Inf) %>%
    rownames_to_column("Symbol") %>%
    select(Symbol, !!paste0("logFC_", label) := logFC, !!paste0("p_", label) := adj.P.Val)
  return(res)
}

# --- 依次计算三个数据集的差异 ---
# 1. RNA-seq
res_rna <- get_degs("Human__CPTAC_COAD__UNC__RNAseq__HiSeq_RNA__03_01_2017__BCM__Gene__BCM_RSEM_UpperQuartile_log2.cct.gz", meta_rna, "RNA")

# 2. Proteome LF (VU平台)
res_lf <- get_degs("Human__CPTAC_COAD__VU__Proteome__QExact__03_01_2017__BCM__Gene__VU_Tumor_LF_UnsharedCounts.cct", meta_rna, "LF")

# 3. Proteome TMT (PNNL平台)
res_tmt <- get_degs("Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct", meta_rna, "TMT")

# 合并三者共同拥有的基因
multi_omice <- res_rna %>%
  inner_join(res_lf, by = "Symbol") %>%
  inner_join(res_tmt, by = "Symbol")

# 计算相关性
cor_rna_lf <- cor(multi_omice$logFC_RNA, multi_omice$logFC_LF, use = "complete.obs")
cor_rna_tmt <- cor(multi_omice$logFC_RNA, multi_omice$logFC_TMT, use = "complete.obs")
cor_lf_tmt <- cor(multi_omice$logFC_LF, multi_omice$logFC_TMT, use = "complete.obs")

library(ggpubr)
library(patchwork)

# 绘图函数
plot_comp <- function(df, x_col, y_col, title) {
  ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(alpha = 0.3, color = "#4682B4") +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
    geom_hline(yintercept = c(-1, 1), linetype = "dotted") +
    stat_cor(method = "pearson") +
    theme_bw() +
    labs(title = title, x = paste("log2FC", x_col), y = paste("log2FC", y_col))
}

p1 <- plot_comp(multi_omice, "logFC_RNA", "logFC_LF", "RNA vs Protein (LF)")
p2 <- plot_comp(multi_omice, "logFC_RNA", "logFC_TMT", "RNA vs Protein (TMT)")
p3 <- plot_comp(multi_omice, "logFC_LF", "logFC_TMT", "Protein (LF) vs Protein (TMT)")

(p1 | p2) / p3 + plot_annotation(title = "Cross-Platform & Cross-Omics Consistency (EMT vs Epi)")

#####其他平台vs#####
library(tidyverse)
library(ggpubr)
library(patchwork)

# 封装一个蛋白组差异分析函数
process_proteome <- function(file_path, meta_data) {
  # 1. 读取
  df <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # 2. 对齐样本
  common <- intersect(colnames(df), meta_data$sample_id)
  df_sub <- df[, common]
  m_sub <- meta_data[match(common, meta_data$sample_id), ]
  # 3. 线性模型
  grp <- factor(m_sub$phenotype, levels = c("Epithelial", "EMT", "Hypermutated"))
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  
  fit <- lmFit(df_sub, design)
  contrast <- makeContrasts(
    EMT_vs_Epi = EMT - Epithelial,
    EMT_vs_Hyper = EMT - Hypermutated,
    Hyper_vs_Epi = Hypermutated - Epithelial,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  return(fit2)
}

# --- 执行蛋白组拟合 (解决你报错找不到对象的问题) ---
lf_fit2 <- process_proteome("Human__CPTAC_COAD__VU__Proteome__QExact__03_01_2017__BCM__Gene__VU_Tumor_LF_UnsharedCounts.cct", meta_sub)
tmt_fit2 <- process_proteome("Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct", meta_sub)

#####绘图#####
# --- 1. 核心提取与绘图函数 ---
plot_multi_omics_panel <- function(rna_f, lf_f, tmt_f, coef_name, title_text) {
  # 提取数据
  get_df <- function(f, lab) {
    topTable(f, coef = coef_name, number = Inf) %>% 
      rownames_to_column("Symbol") %>% 
      select(Symbol, !!paste0("logFC_", lab) := logFC)
  }
  
  res <- get_df(rna_f, "RNA") %>%
    inner_join(get_df(lf_f, "LF"), by = "Symbol") %>%
    inner_join(get_df(tmt_f, "TMT"), by = "Symbol")
  
  # 绘制子图 A: RNA vs LF
  p1 <- ggplot(res, aes(x = logFC_RNA, y = logFC_LF)) +
    geom_point(alpha = 0.3, color = "#4DBFBC", size = 0.8) +
    geom_smooth(method = "lm", color = "darkred", size = 0.5) +
    stat_cor(method = "pearson", size = 3) +
    theme_bw() + labs(subtitle = "RNA vs LF", x = "logFC (RNA)", y = "logFC (LF)")
  
  # 绘制子图 B: RNA vs TMT
  p2 <- ggplot(res, aes(x = logFC_RNA, y = logFC_TMT)) +
    geom_point(alpha = 0.3, color = "#E88471", size = 0.8) +
    geom_smooth(method = "lm", color = "darkred", size = 0.5) +
    stat_cor(method = "pearson", size = 3) +
    theme_bw() + labs(subtitle = "RNA vs TMT", x = "logFC (RNA)", y = "logFC (TMT)")
  
  return((p1 + p2) + plot_annotation(title = title_text))
}

# --- 2. 生成三组对比的面板图 ---

# 第一组: EMT vs Hypermutated
panel_emt_hyper <- plot_multi_omics_panel(fit2, lf_fit2, tmt_fit2, "EMT_vs_Hyper", "Contrast: EMT vs Hypermutated")

# 第二组: Epithelial vs Hypermutated
panel_epi_hyper <- plot_multi_omics_panel(fit2, lf_fit2, tmt_fit2, "Hyper_vs_Epi", "Contrast: Hypermutated vs Epithelial")

# 第三组: EMT vs Epithelial (经典的)
panel_emt_epi <- plot_multi_omics_panel(fit2, lf_fit2, tmt_fit2, "EMT_vs_Epi", "Contrast: EMT vs Epithelial")

# 将三组对比拼接成一张大图
final_all_comparisons <- panel_emt_hyper / panel_epi_hyper / panel_emt_epi

# 打印并保存
print(final_all_comparisons)
ggsave("All_Subtypes_Multi_Omics_Consistency.pdf", final_all_comparisons, width = 10, height = 15)

#######@差异分析 三平台
# 如果没有安装，请运行：
BiocManager::install("GSEABase")
BiocManager::install("GSVA")

library(GSVA)
library(GSEABase)
library(tidyverse)
library(ggpubr)
# 设置环境变量，防止因为警告导致安装中断
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# 核心：在没有任何包被加载的情况下安装
# 这会更新 BiocGenerics 及其依赖项
BiocManager::install(c("BiocGenerics", "GSEABase", "GSVA"), ask = FALSE, update = TRUE)


#### 加载免疫细胞基因集 ####
devtools::install_github("IOBR/IOBR")
#### 1. 免疫浸润估算 (CIBERSORT) ####
library(IOBR)
# 将 log2 数据还原为线性 scale (CIBERSORT 要求)
rna_linear <- 2^rna_sub - 1

# 强制载入 IOBR 内置的 CIBERSORT 签名矩阵
data("lm22", package = "IOBR") 

# 检查 lm22 是否已经出现在环境中了
if(exists("lm22")){
  message("成功找到 lm22 对象！")
} else {
  message("仍然找不到，请检查 IOBR 是否安装正确。")
}
# 使用 deconvo_tme 函数一键调用 CIBERSORT
# method = "cibersort" 内部自动调用 LM22 基因集
tme_cibersort <- deconvo_tme(eset = rna_linear, 
                             method = "cibersort", 
                             arrays = FALSE) # RNA-seq 数据填 FALSE

# 整理数据：合并亚型信息
imm_results <- tme_cibersort %>%
  rename(sample_id = ID) %>%
  inner_join(meta_sub, by = "sample_id")

# 检查一下前几行，确保 22 种细胞数据已生成
head(imm_results)
#### 2. 数据转换与差异分析 ####

# 转换为长格式
imm_long <- imm_results %>%
  pivot_longer(cols = -c(sample_id, phenotype), 
               names_to = "Cell_Type", 
               values_to = "Proportion")

# 绘制差异箱线图
p_imm_diff <- ggplot(imm_long, aes(x = Cell_Type, y = Proportion, fill = phenotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                               "EMT" = "#FC8D62", 
                               "Hypermutated" = "#8DA0CB")) +
  theme_bw() +
  # 计算组间差异 (Kruskal-Wallis 检验)
  stat_compare_means(aes(group = phenotype), label = "p.signif", method = "anova", size = 3) +
  labs(title = "CIBERSORT Immune Cell Proportions across COAD Subtypes",
       x = "", y = "Cell Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "top",
        panel.grid.major.x = element_blank())

print(p_imm_diff)
# 绘制每个样本的堆叠柱状图
p_stack <- imm_long %>%
  ggplot(aes(x = sample_id, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(. ~ phenotype, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(title = "Immune Composition per Sample", x = "Samples", y = "Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7))

print(p_stack)

#####MOFA2（Multi-Omics Factor Analysis v2）
BiocManager::install("MOFA2")
library(MOFA2)
library(tidyverse)

#' 提取并对齐组学数据的通用函数
#' @param file_path 原始数据文件的路径 (.cct 或 .cct.gz)
#' @param meta_data 包含 sample_id 和 phenotype 的临床信息框
#' @param platform_name 平台标签（用于区分日志输出）
#' @return 返回一个清洗好且与 meta 对齐的 matrix

prepare_omics_data <- function(file_path, meta_data, platform_name = "Unknown") {
  
  message(paste0("正在处理平台: ", platform_name, " ..."))
  
  # 1. 读取数据 (自动处理 .gz)
  # row.names = 1 将第一列 (Gene Symbol) 设为行名
  raw_df <- read.table(file_path, header = TRUE, sep = "\t", 
                       row.names = 1, check.names = FALSE, quote = "")
  
  # 2. 统一样本 ID 格式 (可选：部分 CPTAC 数据会将 ID 中的 - 变成 .)
  # colnames(raw_df) <- gsub("\\.", "-", colnames(raw_df))
  
  # 3. 寻找共同样本 (组学矩阵列名 与 临床信息行名的交集)
  common_samples <- intersect(colnames(raw_df), meta_data$sample_id)
  
  if (length(common_samples) == 0) {
    stop(paste0("错误：平台 ", platform_name, " 与临床数据之间没有共同样本！请检查 ID 格式。"))
  }
  
  # 4. 提取子集并排序（确保顺序完全一致）
  # 只取在 meta 中存在的样本，并按 meta 的顺序排列
  raw_df_sub <- raw_df[, common_samples]
  
  # 5. 转换为 Matrix 并处理缺失值 (蛋白组常见 NA)
  # 如果是蛋白组，可能需要去除全为 NA 的行
  raw_matrix <- as.matrix(raw_df_sub)
  raw_matrix <- raw_matrix[rowSums(!is.na(raw_matrix)) > 0, ]
  
  message(paste0("成功！提取了 ", nrow(raw_matrix), " 个特征和 ", ncol(raw_matrix), " 个样本。"))
  
  return(raw_matrix)
}

# 1. 提取 RNA-seq 数据
rna_raw <- prepare_omics_data(
  file_path = "Human__CPTAC_COAD__UNC__RNAseq__HiSeq_RNA__03_01_2017__BCM__Gene__BCM_RSEM_UpperQuartile_log2.cct.gz", 
  meta_data = meta_rna, 
  platform_name = "RNA-seq"
)

# 2. 提取 VU 平台的 LF 蛋白数据
lf_raw <- prepare_omics_data(
  file_path = "Human__CPTAC_COAD__VU__Proteome__QExact__03_01_2017__BCM__Gene__VU_Tumor_LF_UnsharedCounts.cct", 
  meta_data = meta_rna, 
  platform_name = "Proteome_LF"
)

# 3. 提取 PNNL 平台的 TMT 蛋白数据
tmt_raw <- prepare_omics_data(
  file_path = "Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct", 
  meta_data = meta_rna, 
  platform_name = "Proteome_TMT"
)

# 1. 提取共同样本（确保三个平台和临床数据对齐）
common_samples <- intersect(intersect(colnames(rna_raw), colnames(lf_raw)), colnames(tmt_raw))
common_samples <- intersect(common_samples, meta_rna$sample_id)

# 2. 预处理函数：筛选前 2000 个高变基因/蛋白（提高效率）
prepare_layer <- function(df, samples, n_features = 2000) {
  df_sub <- as.matrix(df[, samples])
  # 筛选高变异特征
  vars <- apply(df_sub, 1, var)
  top_features <- names(sort(vars, decreasing = TRUE))[1:min(n_features, nrow(df_sub))]
  return(df_sub[top_features, ])
}

# 构建 MOFA 输入列表
data_list <- list(
  RNA = prepare_layer(rna_raw, common_samples),
  LF  = prepare_layer(lf_raw, common_samples),
  TMT = prepare_layer(tmt_raw, common_samples)
)

# 创建 MOFA 对象
MOFAobj <- create_mofa(data_list)

# 可视化数据概览（查看缺失值和数据分布）
plot_data_overview(MOFAobj)

# 配置训练选项
model_opts <- get_default_model_options(MOFAobj)
model_opts$num_factors <- 15  # 初始设定 15 个因子，模型会自动丢弃无意义的

train_opts <- get_default_training_options(MOFAobj)
train_opts$convergence_mode <- "medium" # 兼顾速度与准确度

# 准备并运行模型
MOFAobj <- prepare_mofa(MOFAobj, model_options = model_opts, training_options = train_opts)
# 在 run_mofa 这一步增加 use_basilisk = TRUE 参数
MOFAobj <- run_mofa(MOFAobj, 
                    outfile = "COAD_MOFA_Model.hdf5", 
                    use_basilisk = TRUE)
# 观察各因子在不同层中的方差解释比例
plot_variance_explained(MOFAobj, x = "view", y = "factor")
# 1. 修改列名以符合 MOFA2 的硬性要求
sample_metadata <- meta_rna %>% 
  filter(sample_id %in% common_samples) %>%
  rename(sample = sample_id) # 将 sample_id 改为 sample

# 2. 注入元数据 (这次不会报错了)
samples_metadata(MOFAobj) <- sample_metadata

# 3. 检查是否注入成功
head(samples_metadata(MOFAobj))
# 观察各因子在 RNA, LF, TMT 中的解释能力
plot_variance_explained(MOFAobj, x = "view", y = "factor")

# 将临床亚型信息加入模型
samples_metadata(MOFAobj) <- sample_metadata

# 绘制因子在不同亚型中的分布（如 Factor 1）
plot_factor(MOFAobj, 
            factors = 1:3, 
            color_by = "phenotype", 
            add_violin = TRUE, 
            dodge = TRUE) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"))
# 查看 Factor 1 中权重最高的特征
plot_factor1<-plot_weights(MOFAobj, view = "RNA", factor = 1, nfeatures = 10)
plot_factor1_line<- plot_top_weights(MOFAobj, view = "TMT", factor = 1, nfeatures = 10)
plot_factor2<-plot_weights(MOFAobj, view = "RNA", factor = 2, nfeatures = 10)
plot_factor2_line<- plot_top_weights(MOFAobj, view = "TMT", factor = 2, nfeatures = 10)

# 绘制 Factor 1 和 Factor 2 的散点图，观察样本聚类情况
plot_factor(MOFAobj, 
            factors = 1:2, 
            color_by = "phenotype",
            shape_by = "phenotype",
            dot_size = 3) +
  scale_color_manual(values = c("Epithelial" = "#66C2A5", 
                                "EMT" = "#FC8D62", 
                                "Hypermutated" = "#8DA0CB"))
###Factor1因子分析
library(MOFA2)
library(tidyverse)

# --- 1. 提取权重表格 ---
# weights 包含了每个基因对 Factor 1 的贡献值（Loading）
weights <- get_weights(MOFAobj, views = "RNA", factors = 1, as.data.frame = TRUE)

# --- 2. 筛选权重绝对值最大的前 50 个基因 ---
# 权重绝对值越大，说明该基因对该因子的驱动作用越强
top50_genes_df <- weights %>%
  mutate(abs_value = abs(value)) %>%
  arrange(desc(abs_value)) %>%
  head(50) %>%
  select(feature, value, abs_value) %>%
  rename(Symbol = feature, Loading = value)

# 打印表格查看
print(top50_genes_df)

# 保存为 CSV 文件
write.csv(top50_genes_df, "Factor1_Top50_Loadings.csv", row.names = FALSE)
# 加载富集分析核心包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# --- 1. 基因 ID 转换 (Symbol 转 ENTREZID) ---
gene_list <- top50_genes_df$Symbol
gene_convert <- bitr(gene_list, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# --- 2. 运行 GO 富集分析 (BP: 生物学过程) ---
ego <- enrichGO(gene          = gene_convert$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# --- 3. 可视化富集结果 ---
# 柱状图
p1 <- barplot(ego, showCategory = 15) + 
  ggtitle("GO Biological Process: Factor 1 Top 50 Genes")

# 气泡图
p2 <- dotplot(ego, showCategory = 15)

print(p1)
print(p2)

# 保存富集结果表格
write.csv(as.data.frame(ego), "Factor1_GO_Enrichment_Results.csv", row.names = FALSE)

###Factor3因子分析
# --- 1. 提取权重表格 ---
# weights 包含了每个基因对 Factor 3 的贡献值（Loading）
weights3 <- get_weights(MOFAobj, views = "RNA", factors = 3, as.data.frame = TRUE)

# --- 2. 筛选权重绝对值最大的前 50 个基因 ---
# 权重绝对值越大，说明该基因对该因子的驱动作用越强
top50_genes_df3 <- weights3 %>%
  mutate(abs_value = abs(value)) %>%
  arrange(desc(abs_value)) %>%
  head(50) %>%
  select(feature, value, abs_value) %>%
  rename(Symbol = feature, Loading = value)

# 打印表格查看
print(top50_genes_df3)

# 保存为 CSV 文件
write.csv(top50_genes_df, "Factor3_Top50_Loadings.csv", row.names = FALSE)
# 加载富集分析核心包

# --- 1. 基因 ID 转换 (Symbol 转 ENTREZID) ---
gene_list3 <- top50_genes_df3$Symbol
gene_convert3 <- bitr(gene_list3, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# --- 2. 运行 GO 富集分析 (BP: 生物学过程) ---
ego3 <- enrichGO(gene          = gene_convert3$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

# --- 3. 可视化富集结果 ---
# 柱状图
p5 <- barplot(ego3, showCategory = 15) + 
  ggtitle("GO Biological Process: Factor 3 Top 50 Genes")

# 气泡图
p7 <- dotplot(ego3, showCategory = 15)

print(p7)
print(p5)

# 保存富集结果表格
write.csv(as.data.frame(ego), "Factor3_GO_Enrichment_Results.csv", row.names = FALSE)
plot_factor(MOFAobj, factors = 1, color_by = "phenotype", add_violin = TRUE)
plot_factor(MOFAobj, factors = 2, color_by = "phenotype", add_violin = TRUE)
# 提取 Factor 分数并进行统计检验
df <- samples_metadata(MOFAobj)
# 假设 Factor1 已经作为一列存在于 metadata 中，如果没有，请先提取
df$Factor1 <- get_factors(MOFAobj, factors = 1)[[1]]
compare_means(Factor1 ~ phenotype, data = df, method = "anova")

# 同时绘制 Factor 2 和 Factor 3 在各亚型中的得分
plot_factor(MOFAobj, 
            factors = c(2, 3), 
            color_by = "phenotype", 
            add_violin = TRUE,
            dodge = TRUE) +
  scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                               "EMT" = "#FC8D62", 
                               "Hypermutated" = "#8DA0CB"))



###Factor2因子分析
# --- 1. 提取权重表格 ---
# weights 包含了每个基因对 Factor 2 的贡献值（Loading）
weights2 <- get_weights(MOFAobj, views = "RNA", factors = 2, as.data.frame = TRUE)

# --- 2. 筛选权重绝对值最大的前 50 个基因 ---
# 权重绝对值越大，说明该基因对该因子的驱动作用越强
top50_genes_df2 <- weights2 %>%
  mutate(abs_value = abs(value)) %>%
  arrange(desc(abs_value)) %>%
  head(50) %>%
  select(feature, value, abs_value) %>%
  rename(Symbol = feature, Loading = value)

# 打印表格查看
print(top50_genes_df2)

# 保存为 CSV 文件
write.csv(top50_genes_df, "Factor2_Top50_Loadings.csv", row.names = FALSE)
# 加载富集分析核心包

# --- 1. 基因 ID 转换 (Symbol 转 ENTREZID) ---
gene_list2 <- top50_genes_df2$Symbol
gene_convert2 <- bitr(gene_list2, fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)

# --- 2. 运行 GO 富集分析 (BP: 生物学过程) ---
ego2 <- enrichGO(gene          = gene_convert2$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)

# --- 3. 可视化富集结果 ---
# 柱状图
p3 <- barplot(ego2, showCategory = 15) + 
  ggtitle("GO Biological Process: Factor 2 Top 50 Genes")

# 气泡图
p4 <- dotplot(ego2, showCategory = 15)

print(p3)
print(p4)

# 保存富集结果表格
write.csv(as.data.frame(ego), "Factor2_GO_Enrichment_Results.csv", row.names = FALSE)
###三平台overlap
# 1. 加载必要的库
library(UpSetR)
library(tidyverse)

# 2. 提取各平台的 Gene ID (行名) 
# 注意：这里直接从你读取的原始矩阵中提取行名
genes_rna <- rownames(rna_raw)
genes_lf  <- rownames(lf_raw)
genes_tmt <- rownames(tmt_raw)

# 3. 统计并打印各平台数量
stats_df <- data.frame(
  Platform = c("RNA-seq", "Proteomics (LF)", "Proteomics (TMT)"),
  Count = c(length(genes_rna), length(genes_lf), length(genes_tmt))
)

cat("各平台检测到的特征数量统计：\n")
print(stats_df)

# 4. 构建 UpSetR 输入列表
upset_list <- list(
  `RNA-seq` = genes_rna,
  `Proteome (LF)` = genes_lf,
  `Proteome (TMT)` = genes_tmt
)

# 5. 绘制 UpSet 图
# nsets: 集合数量; order.by: 按交集频率排序; main.bar.color: 主柱状图颜色
upset(fromList(upset_list), 
      nsets = 3, 
      order.by = "freq", 
      decreasing = TRUE,
      main.bar.color = "#56B4E9", 
      sets.bar.color = "#D55E00",
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.2), # 调整字体大小
      mb.ratio = c(0.6, 0.4)) # 调整上方柱状图与下方点图的比例
############venn图
library(ggvenn)
library(tidyverse)
# 构建列表对象
venn_list <- list(
  `RNA-seq` = genes_rna,
  `Proteome (LF)` = genes_lf,
  `Proteome (TMT)` = genes_tmt
)

# 绘制 Venn 图
# ggvenn 会自动计算交集并绘图，支持设置颜色和透明度
p <- ggvenn(
  venn_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"), # 经典学术配色
  stroke_size = 0.5, 
  set_name_size = 4,
  text_size = 3
) +
  labs(title = "Overlap of Gene Identifiers Across Platforms") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 显示图形
print(p)

# 提取三个平台共有的基因名 (用于后续分析)
common_genes <- intersect(intersect(genes_rna, genes_lf), genes_tmt)
message(paste("三个平台共同检测到的基因数量为:", length(common_genes)))

# 保存共有基因列表
write.table(common_genes, "Common_Genes_Across_3_Platforms.txt", row.names = FALSE, col_names = FALSE, quote = FALSE)