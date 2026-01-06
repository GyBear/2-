# 加载所需包
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(VennDiagram)
library(gridExtra)
library(grid)
library(factoextra)  # 用于PCA可视化
# 设置工作目录（根据你的路径修改）
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")

# ==============================================================================
# 1. 读取所有6个CSV文件
# ==============================================================================
lf_emt_vs_epi     <- read_csv("LF_EMT_vs_Epithelial_Corrected.csv")
lf_hyper_vs_epi   <- read_csv("LF_Hyper_vs_Epithelial_Corrected.csv")
lf_hyper_vs_emt   <- read_csv("LF_Hyper_vs_EMT_Corrected.csv")

tmt_emt_vs_epi    <- read_csv("tmt_EMT_vs_Epithelial_Corrected.csv")
tmt_hyper_vs_epi  <- read_csv("tmt_Hyper_vs_Epithelial_Corrected.csv")
tmt_hyper_vs_emt  <- read_csv("tmt_Hyper_vs_EMT_Corrected.csv")

# 将所有结果放入一个列表，便于循环处理
diff_list <- list(
  LF_EMT_vs_Epi     = lf_emt_vs_epi,
  LF_Hyper_vs_Epi   = lf_hyper_vs_epi,
  LF_Hyper_vs_EMT   = lf_hyper_vs_emt,
  TMT_EMT_vs_Epi    = tmt_emt_vs_epi,
  TMT_Hyper_vs_Epi  = tmt_hyper_vs_epi,
  TMT_Hyper_vs_EMT  = tmt_hyper_vs_emt
)

# ==============================================================================
# 2. 参数设置
# ==============================================================================
sig_threshold <- 0.05
fc_threshold  <- 0.5   # |logFC| > 0.5 为显著（可根据实际情况调整）

# ==============================================================================
# 3. 循环生成火山图 (Volcano Plot)
# ==============================================================================
for(name in names(diff_list)) {
  df <- diff_list[[name]]
  
  # 计算显著蛋白
  sig_prots <- df %>%
    mutate(sig = adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>%
    filter(sig) %>%
    pull(Protein)
  
  common_with_emt_epi <- intersect(sig_prots, diff_list$LF_EMT_vs_Epi %>% 
                                     filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% 
                                     pull(Protein))
  
  top10 <- head(df %>% arrange(adj.P.Val) %>% pull(Protein), 10)
  label_prots <- unique(c(top10, common_with_emt_epi))
  
  volcano <- EnhancedVolcano(df,
                             lab = df$Protein,
                             x = 'logFC',
                             y = 'adj.P.Val',
                             title = paste0(name, ' Differential Expression'),
                             subtitle = paste0('Significant proteins: ', length(sig_prots)),
                             pCutoff = sig_threshold,
                             FCcutoff = fc_threshold,
                             pointSize = 2.5,
                             labSize = 4.0,
                             col = c('grey70', 'forestgreen', 'royalblue', 'red2'),
                             colAlpha = 0.8,
                             cutoffLineType = 'dashed',
                             cutoffLineCol = 'black',
                             selectLab = label_prots,
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             colConnectors = 'black')
  
  # 保存火山图
  ggsave(paste0("Volcano_", gsub(" ", "_", name), ".png"), 
         volcano, width = 12, height = 9, dpi = 300)
  
  print(volcano)  # 在RStudio中显示
}

# ==============================================================================
# 4. 显著蛋白 Venn 图（两平台同一对比的重叠）
# ==============================================================================
# 4.1 EMT vs Epithelial (LF vs TMT)
lf_sig_emt_epi   <- lf_emt_vs_epi %>% filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% pull(Protein)
tmt_sig_emt_epi  <- tmt_emt_vs_epi %>% filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% pull(Protein)

venn_emt_epi <- venn.diagram(
  x = list(LF = lf_sig_emt_epi, TMT = tmt_sig_emt_epi),
  category.names = c("Label-Free", "TMT"),
  filename = NULL,
  main = "Significant Proteins: EMT vs Epithelial",
  lwd = 2, fill = c(alpha("skyblue", 0.4), alpha("salmon", 0.4)),
  cex = 1.5, cat.cex = 1.5
)

png("Venn_EMT_vs_Epithelial.png", width = 800, height = 800)
grid.draw(venn_emt_epi)
dev.off()
grid.draw(venn_emt_epi)

# 4.2 Hypermutated vs Epithelial
lf_sig_hyper_epi   <- lf_hyper_vs_epi %>% filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% pull(Protein)
tmt_sig_hyper_epi  <- tmt_hyper_vs_epi %>% filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% pull(Protein)

venn_hyper_epi <- venn.diagram(
  x = list(LF = lf_sig_hyper_epi, TMT = tmt_sig_hyper_epi),
  category.names = c("Label-Free", "TMT"),
  filename = NULL,
  main = "Significant Proteins: Hypermutated vs Epithelial",
  lwd = 2, fill = c(alpha("skyblue", 0.4), alpha("salmon", 0.4)),
  cex = 1.5, cat.cex = 1.5
)

png("Venn_Hyper_vs_Epithelial.png", width = 800, height = 800)
grid.draw(venn_hyper_epi)
dev.off()
grid.draw(venn_hyper_epi)

# 4.3 Hypermutated vs EMT
lf_sig_hyper_emt   <- lf_hyper_vs_emt %>% filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% pull(Protein)
tmt_sig_hyper_emt  <- tmt_hyper_vs_emt %>% filter(adj.P.Val < sig_threshold & abs(logFC) > fc_threshold) %>% pull(Protein)

venn_hyper_emt <- venn.diagram(
  x = list(LF = lf_sig_hyper_emt, TMT = tmt_sig_hyper_emt),
  category.names = c("Label-Free", "TMT"),
  filename = NULL,
  main = "Significant Proteins: Hypermutated vs EMT",
  lwd = 2, fill = c(alpha("skyblue", 0.4), alpha("salmon", 0.4)),
  cex = 1.5, cat.cex = 1.5
)

png("Venn_Hyper_vs_EMT.png", width = 800, height = 800)
grid.draw(venn_hyper_emt)
dev.off()
grid.draw(venn_hyper_emt)

# ==============================================================================
# 5. PCA 可视化（基于原始预处理后的表达矩阵）
# ==============================================================================
# 注意：PCA 需要原始表达矩阵（样本为行，蛋白为列），而非差异结果
# 这里假设你仍能访问预处理后的表达矩阵 lf_processed$expr 和 tmt_processed$expr
# 以及 pheno（Integrated.Phenotype）
pca_plot_improved <- function(expr_mat, pheno, title_prefix) {
  # expr_mat: 样本为行，蛋白为列（必须确保！）
  # pheno: factor，与行数相同
  
  # 确保 expr_mat 是样本为行
  if(ncol(expr_mat) < nrow(expr_mat)) {
    warning("矩阵可能是转置的，正在自动转置...")
    expr_mat <- t(expr_mat)
  }
  
  # PCA 计算（样本为行）
  pca_res <- prcomp(expr_mat, scale. = TRUE)
  
  # 定义颜色：必须是 named vector，且名字精确匹配 pheno 的 levels
  pheno_levels <- levels(pheno)
  col_vector <- c("Epithelial" = "#1f77b4",    # 深蓝
                  "EMT"        = "#ff7f0e",    # 橙色
                  "Hypermutated" = "#2ca02c")   # 深绿
  
  # 只保留实际出现的 levels 的颜色
  col_vector <- col_vector[names(col_vector) %in% pheno_levels]
  
  p <- fviz_pca_ind(pca_res,
                    geom.ind = "point",
                    col.ind = pheno,                  # factor
                    palette = col_vector,             # named vector，关键！
                    addEllipses = TRUE,
                    ellipse.type = "norm",
                    ellipse.level = 0.85,             # 更大椭圆
                    ellipse.alpha = 0.15,
                    legend.title = "Phenotype",
                    title = paste(title_prefix, "PCA - Samples"),
                    pointsize = 3.5,
                    pointshape = 19,
                    mean.point = FALSE) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  
  return(p)
}

# ==============================================================================
# 生成并保存（确保 lf_processed$expr 和 tmt_processed$expr 存在）
# ==============================================================================
# Label-Free
p_lf_improved <- pca_plot_improved(lf_processed$expr, 
                                   lf_processed$pheno, 
                                   "Label-Free")

print(p_lf_improved)
ggsave("PCA_LabelFree_Improved.png", p_lf_improved, width = 11, height = 8, dpi = 300)

# TMT
p_tmt_improved <- pca_plot_improved(tmt_processed$expr, 
                                    tmt_processed$pheno, 
                                    "TMT")

print(p_tmt_improved)
ggsave("PCA_TMT_Improved.png", p_tmt_improved, width = 11, height = 8, dpi = 300)
# ==============================================================================
# 完成提示
# ==============================================================================
message("所有火山图、Venn图和PCA图已生成并保存到工作目录！")