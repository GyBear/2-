# 加载所需包
library(tidyverse)
library(ggplot2)
library(ggpubr)     
library(EnhancedVolcano)  
library(gridExtra)

# 设置工作目录（根据你的路径修改）
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")

# ==============================================================================
# 1. 读取所有6个差异结果CSV
# ==============================================================================
lf_emt_epi   <- read_csv("LF_EMT_vs_Epithelial_Corrected.csv")
lf_hyper_emt <- read_csv("LF_Hyper_vs_EMT_Corrected.csv")
lf_hyper_epi <- read_csv("LF_Hyper_vs_Epithelial_Corrected.csv")

tmt_emt_epi   <- read_csv("tmt_EMT_vs_Epithelial_Corrected.csv")
tmt_hyper_emt <- read_csv("tmt_Hyper_vs_EMT_Corrected.csv")
tmt_hyper_epi <- read_csv("tmt_Hyper_vs_Epithelial_Corrected.csv")

# ==============================================================================
# 2. 定义显著蛋白阈值
# ==============================================================================
sig_p   <- 0.05      # adj.P.Val < 0.05
sig_fc  <- 0.5       # |logFC| > 0.5

# ==============================================================================
# 3. 超几何检验（Fisher's Exact Test）—— 三组对比的重叠显著性
# ==============================================================================
# 背景：所有在两个平台都检测到的蛋白（取并集）
all_proteins_lf  <- unique(c(lf_emt_epi$Protein, lf_hyper_emt$Protein, lf_hyper_epi$Protein))
all_proteins_tmt <- unique(c(tmt_emt_epi$Protein, tmt_hyper_emt$Protein, tmt_hyper_epi$Protein))
background <- intersect(all_proteins_lf, all_proteins_tmt)  # 共同检测到的蛋白作为背景
N <- length(background)

hypergeom_test <- function(sig_lf, sig_tmt, background) {
  # sig_lf / sig_tmt: 显著蛋白向量
  k <- length(intersect(sig_lf, sig_tmt))           # 重叠显著蛋白数
  n <- length(sig_lf)                               # LF 显著蛋白数
  m <- length(sig_tmt)                              # TMT 显著蛋白数
  K <- length(intersect(sig_lf, background))        # 背景中LF显著蛋白数（通常等于n）
  
  # Fisher's exact test (2x2 列联表)
  contingency <- matrix(c(k, n - k,
                          m - k, N - n - m + k), nrow = 2)
  test <- fisher.test(contingency, alternative = "greater")
  
  list(overlap = k,
       lf_sig = n,
       tmt_sig = m,
       background = N,
       p_value = test$p.value,
       odds_ratio = test$estimate)
}

# 3.1 EMT vs Epithelial
sig_lf_emt_epi  <- lf_emt_epi   %>% filter(adj.P.Val < sig_p, abs(logFC) > sig_fc) %>% pull(Protein)
sig_tmt_emt_epi <- tmt_emt_epi  %>% filter(adj.P.Val < sig_p, abs(logFC) > sig_fc) %>% pull(Protein)
res_emt_epi <- hypergeom_test(sig_lf_emt_epi, sig_tmt_emt_epi, background)
cat("EMT vs Epithelial 重叠显著性:\n")
print(res_emt_epi)

# 3.2 Hypermutated vs EMT
sig_lf_hyper_emt  <- lf_hyper_emt %>% filter(adj.P.Val < sig_p, abs(logFC) > sig_fc) %>% pull(Protein)
sig_tmt_hyper_emt <- tmt_hyper_emt %>% filter(adj.P.Val < sig_p, abs(logFC) > sig_fc) %>% pull(Protein)
res_hyper_emt <- hypergeom_test(sig_lf_hyper_emt, sig_tmt_hyper_emt, background)
cat("\nHypermutated vs EMT 重叠显著性:\n")
print(res_hyper_emt)

# 3.3 Hypermutated vs Epithelial
sig_lf_hyper_epi  <- lf_hyper_epi %>% filter(adj.P.Val < sig_p, abs(logFC) > sig_fc) %>% pull(Protein)
sig_tmt_hyper_epi <- tmt_hyper_epi %>% filter(adj.P.Val < sig_p, abs(logFC) > sig_fc) %>% pull(Protein)
res_hyper_epi <- hypergeom_test(sig_lf_hyper_epi, sig_tmt_hyper_epi, background)
cat("\nHypermutated vs Epithelial 重叠显著性:\n")
print(res_hyper_epi)

# ==============================================================================
# 4. 相关性分析 & 散点图（logFC 一致性）
# ==============================================================================
# ==============================================================================
# 修正后的相关性分析函数
# ==============================================================================
plot_fc_cor_fixed <- function(lf_df, tmt_df, contrast_name) {
  merged <- lf_df %>%
    select(Protein, logFC_LF = logFC) %>%
    inner_join(tmt_df %>% select(Protein, logFC_TMT = logFC), by = "Protein")
  
  n_common <- nrow(merged)
  cat("\n", contrast_name, " 共同蛋白数:", n_common, "\n")
  
  # 计算相关系数（Spearman 警告可忽略）
  pearson_res  <- cor.test(merged$logFC_LF, merged$logFC_TMT, method = "pearson")
  spearman_res <- cor.test(merged$logFC_LF, merged$logFC_TMT, method = "spearman")
  
  cat("Pearson  r =", round(pearson_res$estimate, 3), 
      "p =", format.pval(pearson_res$p.value), "\n")
  cat("Spearman rho =", round(spearman_res$estimate, 3), 
      "p =", format.pval(spearman_res$p.value), "\n")
  
  # 散点图 - 关键修复：避免使用 Inf/-Inf 作为 label 位置
  p <- ggplot(merged, aes(x = logFC_TMT, y = logFC_LF)) +
    geom_point(alpha = 0.6, color = "steelblue", size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    
    # 修正：使用固定位置（右下角），并手动拼接标签文本
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
             label = paste0("Pearson r = ", round(pearson_res$estimate, 3),
                            "\np = ", format.pval(pearson_res$p.value)),
             size = 5, color = "darkred", fontface = "bold") +
    
    labs(title = paste(contrast_name, "- log2FC Correlation (LF vs TMT)"),
         subtitle = paste0("n = ", n_common, " common proteins"),
         x = "TMT log2 Fold Change",
         y = "Label-Free log2 Fold Change") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none")
  
  # 保存
  filename <- paste0("FC_Correlation_", gsub(" vs ", "_vs_", contrast_name), ".png")
  ggsave(filename, p, width = 9, height = 8, dpi = 300)
  
  print(p)
  
  return(merged)
}

# ==============================================================================
# 重新执行三组对比
# ==============================================================================
cor_emt_epi   <- plot_fc_cor_fixed(lf_emt_epi,   tmt_emt_epi,   "EMT vs Epithelial")
cor_hyper_emt <- plot_fc_cor_fixed(lf_hyper_emt, tmt_hyper_emt, "Hypermutated vs EMT")
cor_hyper_epi <- plot_fc_cor_fixed(lf_hyper_epi, tmt_hyper_epi, "Hypermutated vs Epithelial")

# ==============================================================================
# 5. 显著蛋白的 logFC 相关性热图
# ==============================================================================
# 加载所需包
library(tidyverse)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(grid)
library(gridExtra)

# ==============================================================================
# 最终美化版热图函数：解决 Spearman 文字遮盖问题 + 更专业布局
# ==============================================================================
pretty_triangle_heatmap_final <- function(lf_df, tmt_df, contrast_name, 
                                          sig_threshold = 0.05, fc_threshold = 0.5) {
  # 步骤1: 提取共同显著蛋白
  sig_lf  <- lf_df  %>% filter(adj.P.Val < sig_threshold, abs(logFC) > fc_threshold) %>% pull(Protein)
  sig_tmt <- tmt_df %>% filter(adj.P.Val < sig_threshold, abs(logFC) > fc_threshold) %>% pull(Protein)
  common_sig <- intersect(sig_lf, sig_tmt)
  
  if(length(common_sig) == 0) {
    message("No common significant proteins for ", contrast_name)
    return(NULL)
  }
  
  cat("\n", contrast_name, ": 共同显著蛋白数 =", length(common_sig), "\n")
  
  # 步骤2: 构建矩阵（行=蛋白，列=平台）
  merged <- lf_df %>%
    select(Protein, logFC_LF = logFC) %>%
    inner_join(tmt_df %>% select(Protein, logFC_TMT = logFC), by = "Protein") %>%
    filter(Protein %in% common_sig)
  
  mat <- merged %>% column_to_rownames("Protein") %>% as.matrix()
  mat_z <- t(scale(t(mat)))  # 行 Z-score
  
  # 步骤3: 计算 Spearman（基于所有共同蛋白，而非仅显著）
  all_merged <- lf_df %>%
    select(Protein, logFC_LF = logFC) %>%
    inner_join(tmt_df %>% select(Protein, logFC_TMT = logFC), by = "Protein")
  
  spearman_test <- cor.test(all_merged$logFC_LF, all_merged$logFC_TMT, method = "spearman")
  spearman_rho <- round(spearman_test$estimate, 3)
  p_label <- ifelse(spearman_test$p.value < 2.2e-16, "p < 2.2e-16", 
                    paste("p =", format.pval(spearman_test$p.value)))
  
  # 步骤4: 热图主体（RdBu 经典双色）
  my_col <- rev(brewer.pal(11, "RdBu"))
  
  p_heat <- pheatmap(mat_z,
                     color = my_col,
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     show_rownames = length(common_sig) <= 80,   # 太多蛋白时隐藏名字
                     show_colnames = TRUE,
                     fontsize_row = 8,
                     fontsize_col = 14,
                     cellwidth = 28,
                     cellheight = 9,
                     border_color = "white",         # 白色边框，更清爽
                     legend = TRUE,
                     main = paste0(contrast_name, "\nCommon Significant Proteins (Row Z-score logFC)"),
                     silent = TRUE)
  
  # 步骤5: Spearman 文字面板（独立放在热图右上角，避免任何遮盖）
  text_df <- data.frame(
    label = paste0("Spearman ρ = ", spearman_rho, "\n", p_label),
    x = 1, y = 1
  )
  
  p_text <- ggplot(text_df, aes(x = x, y = y, label = label)) +
    geom_text(size = 8, fontface = "bold", color = "darkred", hjust = 0.5, vjust = 0.5) +
    xlim(0.5, 1.5) + ylim(0.5, 1.5) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(30, 30, 30, 30))
  
  # 步骤6: 组合布局（热图左侧 + Spearman 文字右侧上角）
  # 使用 grid.arrange 调整宽度比例，让文字面板更小
  final_plot <- grid.arrange(
    p_heat$gtable,
    p_text,
    ncol = 2,
    widths = c(5, 1),        # 热图占大头，文字占小头
    heights = c(1)
  )
  
  # 添加整体标题（可选）
  title_grob <- textGrob(paste0(contrast_name, " - Consistency between LF and TMT"),
                         gp = gpar(fontsize = 18, fontface = "bold"))
  
  final_with_title <- arrangeGrob(title_grob, final_plot, heights = c(0.1, 0.9))
  
  # 保存高清图
  filename <- paste0("Final_Heatmap_CommonSig_", gsub(" vs ", "_vs_", contrast_name), ".png")
  ggsave(filename, final_with_title, width = 14, height = max(9, length(common_sig)/7 + 3), 
         dpi = 400, bg = "white")
  
  grid.draw(final_with_title)
  
  return(final_with_title)
}

# ==============================================================================
# 对三组对比运行最终美化版热图
# ==============================================================================
pretty_triangle_heatmap_final(lf_emt_epi,   tmt_emt_epi,   "EMT vs Epithelial")
pretty_triangle_heatmap_final(lf_hyper_emt, tmt_hyper_emt, "Hypermutated vs EMT")
pretty_triangle_heatmap_final(lf_hyper_epi, tmt_hyper_epi, "Hypermutated vs Epithelial")

# ==============================================================================
# 完成
# ==============================================================================
message("最终美化热图已生成！")
message("改进点：")
message("- Spearman ρ 值移至独立右侧面板，完全避免遮盖")
message("- 白色边框 + 更大间距，更清爽专业")
message("- 整体标题居中，适合论文发表")
message("- 自适应高度，蛋白多时自动拉长")

# ==============================================================================
# 完成
# ==============================================================================
message("所有验证性统计完成！")
message("超几何检验结果已打印，相关性散点图已保存为 FC_Correlation_*.png")