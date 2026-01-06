library(tidyverse)
library(ggpubr)

# 1. 加载与清洗 (保持之前逻辑，增加基因因子排序以固定位置)
target_genes <- c("STAT1", "LTF", "CD14", "C3", "GBP1", "CAMP")

process_data <- function(df, platform) {
  df %>%
    select(sample_id, phenotype, any_of(target_genes)) %>%
    pivot_longer(cols = -c(sample_id, phenotype), names_to = "Gene", values_to = "Expression") %>%
    mutate(Platform = platform)
}

combined_data <- bind_rows(process_data(lf_df, "LF"), process_data(tmt_df, "TMT")) %>%
  filter(!is.na(Expression)) %>%
  mutate(
    phenotype = factor(phenotype, levels = c("Epithelial", "EMT", "Hypermutated")),
    Gene = factor(Gene, levels = target_genes),
    # 创建一个联合标签，方便排列
    Facet_Label = paste0(Gene, " (", Platform, ")")
  )

# 2. 定义比较对
my_comparisons <- list(c("Epithelial", "EMT"), c("EMT", "Hypermutated"), c("Epithelial", "Hypermutated"))

# 3. 绘图：采用 facet_wrap 并在两列分布，避免垂直过度压缩
p <- ggplot(combined_data, aes(x = phenotype, y = Expression, fill = phenotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, size = 0.3) +
  geom_jitter(aes(color = phenotype), width = 0.15, size = 0.2, alpha = 0.4) +
  # 关键修改：使用 facet_wrap 设为 3 列或 4 列，或保持 grid 但大幅增加高度
  # 这里建议按基因分面，平台作为次级分类
  facet_wrap(~ Gene + Platform, scales = "free_y", ncol = 4) + 
  theme_bw() +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  # 显著性标记优化：减小间距
  stat_compare_means(comparisons = my_comparisons, 
                     label = "p.signif", 
                     method = "wilcox.test",
                     step.increase = 0.08, # 减小横线间的间距
                     tip.length = 0.01,
                     size = 3) +
  # 整体 P 值移至左上角，避免与横线冲突
  stat_compare_means(label.y.npc = "top", label.x.npc = 0.05, size = 2.5) +
  labs(x = NULL, y = "Intensity (log-scale)") +
  theme(
    strip.background = element_rect(fill = "white", color = "grey80"),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.y = element_text(size = 10),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"), # 增加子图间距
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# 4. 保存：对于 12 个子图，建议宽度加宽，高度适中
# 若使用 ncol=4 (3行4列)，建议宽 12 高 9
ggsave("Immune_Genes_Optimized.pdf", width = 12, height = 9, device = "pdf", limitsize = FALSE)

library(tidyverse)
library(ggpubr)

# 1. 加载 RNA-seq 数据和临床分型数据
# 注意：RNA-seq 数据通常是基因在行，样本在列，需要转置
rna_raw <- read.table("Human__CPTAC_COAD__UNC__RNAseq__HiSeq_RNA__03_01_2017__BCM__Gene__BCM_RSEM_UpperQuartile_log2.cct.gz", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
clinical <- read.table("Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi", 
                       header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 2. 提取分型信息并清洗
phenotype_info <- as.data.frame(t(clinical["Integrated.Phenotype", ])) %>%
  rownames_to_column("sample_id") %>%
  rename(phenotype = Integrated.Phenotype) %>%
  filter(!is.na(phenotype))

# 3. 定义目标基因并提取数据
# 包含：粘附 (FN1, TLN1, TNC), 复制 (MCM2, MCM4, PCNA, FEN1), 免疫 (CD4, CD14, STAT1)
# 注：MCM 是一个家族，这里选取常用的 MCM2 和 MCM4
target_genes <- c("FN1", "TLN1", "TNC", "MCM2", "MCM4", "PCNA", "FEN1", "CD4", "CD14", "STAT1")

# 提取存在于数据中的基因
present_genes <- intersect(target_genes, rownames(rna_raw))

rna_data <- as.data.frame(t(rna_raw[present_genes, ])) %>%
  rownames_to_column("sample_id")

# 4. 合并数据并添加通路分类标签
plot_data <- inner_join(phenotype_info, rna_data, by = "sample_id") %>%
  pivot_longer(cols = all_of(present_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(
    phenotype = factor(phenotype, levels = c("Epithelial", "EMT", "Hypermutated")),
    # 为基因添加功能分组标签
    Category = case_when(
      Gene %in% c("FN1", "TLN1", "TNC") ~ "Cell-Matrix Adhesion",
      Gene %in% c("MCM2", "MCM4", "PCNA", "FEN1") ~ "DNA Replication",
      Gene %in% c("CD4", "CD14", "STAT1") ~ "Immune Infiltration"
    ),
    Category = factor(Category, levels = c("Cell-Matrix Adhesion", "DNA Replication", "Immune Infiltration"))
  )

# 5. 定义比较组
my_comparisons <- list(c("Epithelial", "EMT"), c("EMT", "Hypermutated"), c("Epithelial", "Hypermutated"))

# 6. 绘图
p <- ggplot(plot_data, aes(x = phenotype, y = Expression, fill = phenotype)) +
  # 绘制半透明箱线图
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, color = "black", size = 0.4) +
  # 添加小抖动点
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.3, aes(color = phenotype)) +
  # 分面：按基因分面，允许 y 轴独立缩放
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  theme_pubr() +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  # 添加两两比较的显著性 (星号表示)
  stat_compare_means(comparisons = my_comparisons, 
                     label = "p.signif", 
                     method = "wilcox.test",
                     step.increase = 0.1,
                     tip.length = 0.01) +
  # 添加总体的 Kruskal-Wallis 检验
  stat_compare_means(label.y.npc = "top", label.x.npc = 0.05, size = 3) +
  labs(title = "RNA-seq Expression: Adhesion, Replication & Immune Genes",
       x = NULL, 
       y = "Log2 (UpperQuartile Normalized)") +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold.italic", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "none"
  )

# 7. 保存结果
print(p)
ggsave("RNAseq_Immune_Adhesion_Replication.pdf", width = 12, height = 10)








