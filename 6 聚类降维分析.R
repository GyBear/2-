library(tidyverse)
library(ggplot2)
library(ggfortify) # 简化 PCA 绘图的包
library(patchwork) # 用于合并 TMT 和 LF 的图片

# --- 1. 处理临床元数据
raw_clin <- read.delim("Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi", 
                       header = FALSE)
t_clin <- t(raw_clin)
colnames(t_clin) <- t_clin[1, ]
clin_df <- as.data.frame(t_clin[-1, ], stringsAsFactors = FALSE)

sample_metadata <- clin_df %>%
  dplyr::select(sample_id = attrib_name, phenotype = Integrated.Phenotype) %>%
  filter(phenotype %in% c("Epithelial", "EMT", "Hypermutated"))
plot_custom_pca <- function(expr_matrix, meta_df, platform_name, is_log_scale = FALSE) {
  
  # 1. 样本对齐
  common_samples <- intersect(colnames(expr_matrix), meta_df$sample_id)
  dat <- expr_matrix[, common_samples]
  meta <- meta_df[match(common_samples, meta_df$sample_id), ]
  
  # 2. 预处理
  # 过滤缺失值过多的蛋白 (PCA 对 NA 非常敏感)
  dat <- dat[rowSums(is.na(dat)) < ncol(dat) * 0.2, ] 
  
  # 均值填充 NA
  dat_imputed <- t(apply(dat, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  
  # LF 数据需要 Log2 转换
  if(!is_log_scale) dat_imputed <- log2(dat_imputed + 1)
  
  # 3. 执行 PCA (prcomp 默认对行进行主成分提取，所以要转置，让行为样本)
  pca_res <- prcomp(t(dat_imputed), scale. = TRUE)
  
  # 4. 绘图
  # 使用 autoplot 自动提取 PC1, PC2 并结合元数据
  p <- autoplot(pca_res, data = meta, colour = 'phenotype', frame = TRUE, frame.type = 'norm') +
    scale_color_manual(values = c("Epithelial" = "#66C2A5", 
                                  "EMT" = "#FC8D62", 
                                  "Hypermutated" = "#8DA0CB")) +
    scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                                 "EMT" = "#FC8D62", 
                                 "Hypermutated" = "#8DA0CB")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = paste0(platform_name, " Global PCA"),
         subtitle = "Clustering of CRC Subtypes",
         colour = "Subtype", fill = "Subtype")
  
  return(p)
}
# 读取数据
tmt_raw <- read.table("Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
lf_raw <- read.table("Human__CPTAC_COAD__VU__Proteome__QExact__03_01_2017__BCM__Gene__VU_Tumor_LF_UnsharedCounts.cct", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 绘图
pca_tmt <- plot_custom_pca(tmt_raw, sample_metadata, "TMT", is_log_scale = TRUE)
pca_lf <- plot_custom_pca(lf_raw, sample_metadata, "LF", is_log_scale = FALSE)

# 合并出图
final_plot <- pca_tmt + pca_lf + plot_layout(guides = 'collect')
ggsave("Figure_Global_PCA_Subtypes.pdf", final_plot, width = 12, height = 5)

print(final_plot)

####T-SNE analysis###
library(Rtsne)     # t-SNE 核心包
library(tidyverse)
library(ggplot2)
library(patchwork)

# --- 1. 样本元数据对齐 (确保 clinical 数据已转置处理) ---
# 这里沿用之前处理好的 sample_metadata (包含 sample_id 和 phenotype)
plot_tsne_comparison <- function(expr_matrix, meta_df, platform_name, is_log_scale = FALSE, perplexity = 15) {
  
  # 1. 样本与元数据对齐
  common_samples <- intersect(colnames(expr_matrix), meta_df$sample_id)
  dat <- expr_matrix[, common_samples]
  meta <- meta_df[match(common_samples, meta_df$sample_id), ]
  
  # 2. 预处理
  # 剔除缺失值过多的行 (t-SNE 不允许 NA)
  dat <- dat[rowSums(is.na(dat)) < ncol(dat) * 0.2, ]
  
  # 均值填充 NA
  dat_imputed <- t(apply(dat, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  
  # LF 数据 Log 转换
  if(!is_log_scale) dat_imputed <- log2(dat_imputed + 1)
  
  # 3. 运行 t-SNE
  # 注意：Rtsne 输入的是样本为行、蛋白为列的矩阵
  # check_duplicates = FALSE 是因为蛋白质表达谱通常没有完全重复的样本
  set.seed(123) # t-SNE 是随机算法，固定种子保证结果可重复
  tsne_res <- Rtsne(t(dat_imputed), dims = 2, perplexity = perplexity, 
                    verbose = FALSE, check_duplicates = FALSE)
  
  # 4. 整理绘图数据
  plot_dat <- data.frame(
    tsne1 = tsne_res$Y[,1],
    tsne2 = tsne_res$Y[,2],
    Subtype = meta$phenotype
  )
  
  # 5. 绘图
  p <- ggplot(plot_dat, aes(x = tsne1, y = tsne2, color = Subtype, fill = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = Subtype)) + # 添加置信椭圆
    scale_color_manual(values = c("Epithelial" = "#66C2A5", 
                                  "EMT" = "#FC8D62", 
                                  "Hypermutated" = "#8DA0CB")) +
    scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                                 "EMT" = "#FC8D62", 
                                 "Hypermutated" = "#8DA0CB")) +
    theme_minimal() +
    labs(title = paste0(platform_name, " t-SNE Analysis"),
         subtitle = paste0("Perplexity = ", perplexity),
         x = "t-SNE dimension 1", y = "t-SNE dimension 2") +
    theme(panel.grid = element_blank(), axis.line = element_line(color = "black"))
  
  return(p)
}
# 运行 TMT 平台 t-SNE
p_tsne_tmt <- plot_tsne_comparison(tmt_raw, sample_metadata, "TMT", is_log_scale = TRUE, perplexity = 20)

# 运行 LF 平台 t-SNE
p_tsne_lf <- plot_tsne_comparison(lf_raw, sample_metadata, "LF", is_log_scale = FALSE, perplexity = 20)

# 合并展示
tsne_final <- p_tsne_tmt + p_tsne_lf + plot_layout(guides = 'collect')
ggsave("Figure_tSNE_Subtypes_Comparison.pdf", tsne_final, width = 12, height = 5)

print(tsne_final)
####UMAP analysis###
library(umap)
library(tidyverse)
library(ggplot2)
library(patchwork)

# 确保之前的 sample_metadata 已经加载，且列名为 sample_id 和 phenotype
plot_umap_comparison <- function(expr_matrix, meta_df, platform_name, is_log_scale = FALSE) {
  
  # 1. 样本与元数据对齐
  common_samples <- intersect(colnames(expr_matrix), meta_df$sample_id)
  dat <- expr_matrix[, common_samples]
  meta <- meta_df[match(common_samples, meta_df$sample_id), ]
  
  # 2. 预处理 (UMAP 不支持 NA)
  dat <- dat[rowSums(is.na(dat)) < ncol(dat) * 0.2, ]
  dat_imputed <- t(apply(dat, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  
  if(!is_log_scale) dat_imputed <- log2(dat_imputed + 1)
  
  # 3. 运行 UMAP
  # UMAP 同样要求样本为行，蛋白为列
  set.seed(42) # 固定种子保证可重复性
  umap_config <- umap.defaults
  umap_config$n_neighbors <- 15 # 邻居数，可根据样本量微调
  
  umap_res <- umap(t(dat_imputed), config = umap_config)
  
  # 4. 整理绘图数据
  plot_dat <- data.frame(
    UMAP1 = umap_res$layout[,1],
    UMAP2 = umap_res$layout[,2],
    Subtype = meta$phenotype
  )
  
  # 5. 绘图
  p <- ggplot(plot_dat, aes(x = UMAP1, y = UMAP2, color = Subtype, fill = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = Subtype)) +
    scale_color_manual(values = c("Epithelial" = "#66C2A5", 
                                  "EMT" = "#FC8D62", 
                                  "Hypermutated" = "#8DA0CB")) +
    scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                                 "EMT" = "#FC8D62", 
                                 "Hypermutated" = "#8DA0CB")) +
    theme_bw() +
    labs(title = paste0(platform_name, " UMAP Projection"),
         subtitle = "Non-linear manifold visualization",
         x = "UMAP dimension 1", y = "UMAP dimension 2") +
    theme(panel.grid = element_blank(),
          legend.position = "right")
  
  return(p)
}
# 运行 TMT 平台 UMAP
p_umap_tmt <- plot_umap_comparison(tmt_raw, sample_metadata, "TMT", is_log_scale = TRUE)

# 运行 LF 平台 UMAP
p_umap_lf <- plot_umap_comparison(lf_raw, sample_metadata, "LF", is_log_scale = FALSE)

# 合并展示
umap_final <- p_umap_tmt + p_umap_lf + plot_layout(guides = 'collect')
ggsave("Figure_UMAP_Subtypes_Comparison.pdf", umap_final, width = 12, height = 5)

print(umap_final)
####PLS-DA analysis####
BiocManager::install("mixOmics")
library(mixOmics)
library(tidyverse)
library(ggplot2)
library(patchwork)
plot_plsda_comparison <- function(expr_matrix, meta_df, platform_name, is_log_scale = FALSE) {
  
  # 1. 样本与元数据对齐
  common_samples <- intersect(colnames(expr_matrix), meta_df$sample_id)
  dat <- expr_matrix[, common_samples]
  meta <- meta_df[match(common_samples, meta_df$sample_id), ]
  
  # 2. 预处理
  # 过滤缺失值 > 20% 的蛋白
  dat <- dat[rowSums(is.na(dat)) < ncol(dat) * 0.2, ]
  
  # 均值填充 NA (mixOmics 要求矩阵完整)
  dat_imputed <- t(apply(dat, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  
  if(!is_log_scale) dat_imputed <- log2(dat_imputed + 1)
  
  # 3. 运行 PLS-DA
  # X: 样本为行，蛋白为列; Y: 分组标签
  X <- t(dat_imputed)
  Y <- as.factor(meta$phenotype)
  
  plsda_res <- plsda(X, Y, ncomp = 2) # 计算前两个主成分
  
  # 4. 绘图 (提取坐标)
  plot_dat <- as.data.frame(plsda_res$variates$X)
  plot_dat$Subtype <- Y
  
  # 计算解释方差百分比
  expl_var <- round(plsda_res$prop_expl_var$X * 100, 1)
  
  p <- ggplot(plot_dat, aes(x = comp1, y = comp2, color = Subtype, fill = Subtype)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = Subtype)) +
    scale_color_manual(values = c("Epithelial" = "#66C2A5", 
                                  "EMT" = "#FC8D62", 
                                  "Hypermutated" = "#8DA0CB")) +
    scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                                 "EMT" = "#FC8D62", 
                                 "Hypermutated" = "#8DA0CB")) +
    theme_bw() +
    labs(title = paste0(platform_name, " PLS-DA Analysis"),
         x = paste0("Component 1 (", expl_var[1], "%)"),
         y = paste0("Component 2 (", expl_var[2], "%)")) +
    theme(panel.grid = element_blank())
  
  return(list(plot = p, model = plsda_res))
}
# 运行分析
res_tmt <- plot_plsda_comparison(tmt_raw, sample_metadata, "TMT", is_log_scale = TRUE)
res_lf <- plot_plsda_comparison(lf_raw, sample_metadata, "LF", is_log_scale = FALSE)

# 合并图片展示
plsda_final <- res_tmt$plot + res_lf$plot + plot_layout(guides = 'collect')
print(plsda_final)
ggsave("Figure_PLSDA_Subtypes.pdf", plsda_final, width = 12, height = 5)
# 计算 VIP 分数
vip_scores <- vip(res_tmt$model)
# 提取 Component 1 中 VIP 值最高的前 20 个蛋白
top_vip_proteins <- sort(vip_scores[, 1], decreasing = TRUE) %>% head(20)

# 可视化 VIP 蛋白
vip_df <- data.frame(Protein = names(top_vip_proteins), VIP = as.numeric(top_vip_proteins))
ggplot(vip_df, aes(x = reorder(Protein, VIP), y = VIP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 VIP Proteins (TMT Component 1)", x = "Protein Symbol", y = "VIP Score")
