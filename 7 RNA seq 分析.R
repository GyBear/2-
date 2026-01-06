####RNA Seq analysis####
# 设置工作目录（根据你的路径修改）
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")
library(tidyverse)
library(limma)

# --- 1. 读取转录组数据 ---
# 该数据已经是 log2 处理过的 RSEM UpperQuartile 值
rna_raw <- read.table("Human__CPTAC_COAD__UNC__RNAseq__HiSeq_RNA__03_01_2017__BCM__Gene__BCM_RSEM_UpperQuartile_log2.cct.gz", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# --- 2. 处理临床亚型数据 ---
raw_clin <- read.delim("Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi", header = FALSE)
t_clin <- t(raw_clin)
colnames(t_clin) <- t_clin[1, ]
clin_df <- as.data.frame(t_clin[-1, ], stringsAsFactors = FALSE)

# 提取亚型并过滤
meta_rna <- clin_df %>%
  dplyr::select(sample_id = attrib_name, phenotype = Integrated.Phenotype) %>%
  filter(phenotype %in% c("Epithelial", "EMT", "Hypermutated"))

# --- 3. 样本对齐 ---
common_samples <- intersect(colnames(rna_raw), meta_rna$sample_id)
rna_sub <- rna_raw[, common_samples]
meta_sub <- meta_rna[match(common_samples, meta_rna$sample_id), ]

message("对齐完成，共包含 ", length(common_samples), " 个样本。")

# 确保分组为因子且顺序固定
group <- factor(meta_sub$phenotype, levels = c("Epithelial", "EMT", "Hypermutated"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 线性模型拟合
fit <- lmFit(rna_sub, design)

# --- 核心步骤：设置三个两两比较的对比矩阵 ---
contrast.matrix <- makeContrasts(
  EMT_vs_Epi = EMT - Epithelial,
  EMT_vs_Hyper = EMT - Hypermutated,
  Hyper_vs_Epi = Hypermutated - Epithelial,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
####火山图分析####
library(ggrepel)
library(patchwork)

draw_volcano <- function(fit, coef_name, title) {
  # 1. 提取差异分析结果
  res <- topTable(fit, coef = coef_name, number = Inf) %>%
    rownames_to_column("Symbol")
  
  # 2. 标记上下调
  # 阈值设定：|logFC| > 1 且 adj.P.Val < 0.05
  res <- res %>%
    mutate(change = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "UP",
      adj.P.Val < 0.05 & logFC < -1 ~ "DOWN",
      TRUE ~ "NOT"
    ))
  
  # 3. 准备标注最显著的前10个基因
  top10 <- res %>% filter(change != "NOT") %>% arrange(adj.P.Val) %>% head(10)
  
  # 4. 绘图
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_color_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NOT" = "grey")) +
    geom_text_repel(data = top10, aes(label = Symbol), size = 3, max.overlaps = 20) +
    # 添加辅助线
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = "log2 Fold Change", y = "-log10 adj.P-value") +
    theme(panel.grid = element_blank(), legend.position = "none")
  
  return(p)
}
# 分别生成三张图
p1 <- draw_volcano(fit2, "EMT_vs_Epi", "EMT vs Epithelial")
p2 <- draw_volcano(fit2, "EMT_vs_Hyper", "EMT vs Hypermutated")
p3 <- draw_volcano(fit2, "Hyper_vs_Epi", "Hyper vs Epithelial")

# 使用 patchwork 横向拼接
final_volcano <- p1 + p2 + p3 + plot_annotation(
  title = "Differential Gene Expression across COAD Subtypes",
  subtitle = "Cutoff: |log2FC| > 1 & adj.P < 0.05",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

# 打印并保存
print(final_volcano)
ggsave("CRC_Subtypes_Volcano_Plots.pdf", final_volcano, width = 15, height = 5)

####韦恩图分析####
library(ggVennDiagram)

# --- 1. 提取三组对比的显著差异基因名单 (adj.P < 0.05 & |logFC| > 1) ---

get_deg_list <- function(fit, coef_name) {
  res <- topTable(fit, coef = coef_name, number = Inf)
  genes <- rownames(res[res$adj.P.Val < 0.05 & abs(res$logFC) > 1, ])
  return(genes)
}

deg_EMT_vs_Epi <- get_deg_list(fit2, "EMT_vs_Epi")
deg_EMT_vs_Hyper <- get_deg_list(fit2, "EMT_vs_Hyper")
deg_Hyper_vs_Epi <- get_deg_list(fit2, "Hyper_vs_Epi")

# 将它们整合到一个 List 中
venn_list <- list(
  EMT_vs_Epithelial = deg_EMT_vs_Epi,
  EMT_vs_Hypermutated = deg_EMT_vs_Hyper,
  Hyper_vs_Epithelial = deg_Hyper_vs_Epi
)
# --- 2. 绘制韦恩图 ---
p_venn <- ggVennDiagram(venn_list, 
                        label_alpha = 0, # 去除标签背景
                        edge_size = 0.5) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + # 颜色梯度
  theme(legend.position = "none") +
  labs(title = "Overlap of Differentially Expressed Genes",
       subtitle = "RNA-seq Level (|log2FC| > 1 & adj.P < 0.05)")

print(p_venn)
ggsave("CRC_Subtypes_Venn_Diagram.pdf", p_venn, width = 8, height = 7)
# 提取三组对比共同拥有的差异基因 (交集)
common_degs <- Reduce(intersect, venn_list)

message("三组对比共有的核心差异基因数量：", length(common_degs))
print(common_degs)

# 将交集基因保存，方便做下一步的功能富集分析
write.table(common_degs, "Common_Core_DEGs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#### PCA 分析 ####
library(ggfortify) # 方便使用 autoplot 绘制 PCA

# 1. 准备数据：PCA 要求样本为行，基因为列
# 我们使用对齐后的 rna_sub 数据
# 注意：通常 PCA 会选择表达量波动较大的基因（高变基因）来提高分辨率
# 这里我们计算每个基因的方差，选取前 2000 个高变基因
rv <- apply(rna_sub, 1, var)
select_genes <- names(sort(rv, decreasing = TRUE))[1:2000]
pca_data <- t(rna_sub[select_genes, ])

# 2. 执行 PCA
# center = TRUE: 中心化（减去均值）
# scale. = TRUE: 标准化（除以标准差），推荐执行以消除量纲影响
pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# 3. 绘制 PCA 图
p_pca <- autoplot(pca_res, 
                  data = meta_sub, 
                  colour = 'phenotype', 
                  frame = TRUE,           # 添加置信椭圆
                  frame.type = 'norm',    # 假设正态分布
                  size = 3) +
  scale_color_manual(values = c("Epithelial" = "#66C2A5", 
                                "EMT" = "#FC8D62", 
                                "Hypermutated" = "#8DA0CB")) +
  scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                               "EMT" = "#FC8D62", 
                               "Hypermutated" = "#8DA0CB")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "PCA Plot of RNA-seq Data",
       subtitle = "Based on Top 2000 Highly Variable Genes",
       colour = "Subtype",
       fill = "Subtype")

# 4. 展示并保存
print(p_pca)
ggsave("CRC_RNAseq_PCA.pdf", p_pca, width = 8, height = 6)

library(ggpubr)

#### 核心交叉基因表达可视化 (含组间显著性标识) ####

# 1. 定义需要相互比较的分组对
# 这将在图上自动绘制比较线和显著性符号
my_comparisons <- list( 
  c("Epithelial", "EMT"), 
  c("Epithelial", "Hypermutated"), 
  c("EMT", "Hypermutated") 
)

# 2. 提取数据 (保持你原来的逻辑)
plot_genes <- head(common_degs, 9) 
expr_long <- rna_sub[plot_genes, ] %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "sample_id", values_to = "Expression") %>%
  left_join(meta_sub, by = "sample_id")

# 3. 绘图
p_box_sig <- ggplot(expr_long, aes(x = phenotype, y = Expression, fill = phenotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
  facet_wrap(~Gene, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = c("Epithelial" = "#66C2A5", 
                               "EMT" = "#FC8D62", 
                               "Hypermutated" = "#8DA0CB")) +
  theme_bw() +
  # --- 添加全局 ANOVA P值 (放在子图最上方) ---
  stat_compare_means(method = "anova", label.y.npc = "top", size = 3) + 
  # --- 添加两两比较的显著性标识 ---
  # label = "p.signif" 会显示星号 (*, **, ***)
  # 若想显示具体数字 P 值，可改为 label = "p.format"
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     step.increase = 0.1) + # 自动调整连线高度，防止重叠
  labs(title = "Expression Profiles of Core Intersection Genes",
       subtitle = "RNA-seq Level with Pairwise Comparisons",
       x = "Molecular Subtype",
       y = "log2(RSEM UpperQuartile)") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"))

# 4. 展示并保存
print(p_box_sig)
ggsave("Core_Genes_Boxplot_Significance.pdf", p_box_sig, width = 10, height = 12)

####heatmap analysis####
#### 补充差异提取 ####
# 使用 coef = 1:3 或者是默认不指定 coef，topTable 会对 design 中的所有比较进行 F-test
# 这将选出在三组亚型中表现不一致（即有显著差异）的基因
anova_rna <- topTable(fit2, number = Inf, adjust.method = "BH")

# 检查一下前几个基因
head(anova_rna)

library(pheatmap)

#### RNA-seq 热图绘制 ####

# 1. 筛选 P 值最小的前 100 个基因名
top100_genes <- rownames(anova_rna)[1:100]

# 2. 提取表达矩阵并进行 Z-score 标准化
# 为什么要标准化：基因间基础表达量差异巨大，标准化后可观察相对变化趋势
heatmap_mat <- rna_sub[top100_genes, ]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# 3. 准备列注释（样本分组标签）
# 使用 remove_rownames() 清除所有现有的行名索引
annotation_col <- meta_sub %>%
  dplyr::select(sample_id, phenotype) %>%
  remove_rownames() %>%               # 先清空行名
  column_to_rownames("sample_id")     # 再设置新的行名

# 4. 准备颜色配置（与之前图表颜色保持一致）
ann_colors <- list(
  phenotype = c("Epithelial" = "#66C2A5", 
                "EMT" = "#FC8D62", 
                "Hypermutated" = "#8DA0CB")
)

# 5. 绘制热图
# 如果基因名太挤，可以将 show_rownames 设为 FALSE
p_heatmap <- pheatmap(heatmap_mat_scaled,
                      cluster_rows = TRUE,        # 对基因进行聚类
                      cluster_cols = TRUE,        # 对样本进行聚类
                      show_colnames = FALSE,      # 样本名通常较多，建议隐藏
                      show_rownames = TRUE,       # 显示前100个基因名
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors,
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 经典红白蓝配色
                      fontsize_row = 7,           # 调小基因名字号以适应版面
                      main = "Top 100 DEGs across COAD Subtypes (ANOVA)",
                      border_color = NA)          # 去掉格子边框更美观

# 6. 保存热图
# pheatmap 保存有时需要用 pdf() 函数包裹
pdf("CRC_Top100_Heatmap.pdf", width = 10, height = 12)
grid::grid.newpage()
grid::grid.draw(p_heatmap$gtable)
dev.off()

#### 样本匹配严格自检 ####

# 1. 检查数量是否相等
count_match <- ncol(rna_sub) == nrow(meta_sub)

# 2. 检查所有的表达矩阵样本 ID 是否都在临床信息中
all_included <- all(colnames(rna_sub) %in% meta_sub$sample_id)

# 3. 检查【顺序】是否完全一致 (这是热图正确的最核心前提)
identical_order <- all(colnames(rna_sub) == meta_sub$sample_id)

# 输出自检报告
message("--- 样本对齐自检报告 ---")
message("1. 数量一致性: ", count_match)
message("2. ID 完整性: ", all_included)
message("3. 顺序一致性: ", identical_order)

# 4. 如果顺序不一致，执行强制重排序（安全保险）
if (!identical_order) {
  message("警告：顺序不一致！正在重新按矩阵列名对齐临床信息...")
  meta_sub <- meta_sub[match(colnames(rna_sub), meta_sub$sample_id), ]
  # 再次确认
  final_check <- all(colnames(rna_sub) == meta_sub$sample_id)
  message("强制对齐结果: ", final_check)
} else {
  message("恭喜：样本顺序完美匹配，可以直接进行热图绘制。")
}