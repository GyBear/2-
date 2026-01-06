library(tidyverse)
library(ComplexHeatmap) # 绘图核心包
library(circlize)       # 颜色控制
library(impute)         # 处理缺失值
library(dplyr)
# 设置工作目录（根据你的路径修改）
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")
# --- 1. 读取原始表达矩阵 (.cct 文件通常是 Tab 分隔) ---
# TMT 数据 (Log Ratio)
tmt_raw <- read.table("Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# LF 数据 (Counts - 需要 Log 处理)
lf_raw <- read.table("Human__CPTAC_COAD__VU__Proteome__QExact__03_01_2017__BCM__Gene__VU_Tumor_LF_UnsharedCounts.cct", 
                     header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# --- 2. 准备样本元数据 (Metadata) ---
# --- 2.1. 读取并清洗临床元数据 (.tsi) ---
# tsi 文件通常没有标准表头，需要转置处理
raw_clin <- read.delim("Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi", 
                       header = FALSE, stringsAsFactors = FALSE)

# 转置矩阵：行变列，列变行
t_clin <- t(raw_clin)
# 将第一行（原本的属性名）设置为列名
colnames(t_clin) <- t_clin[1, ]
# 去除第一行，转为数据框
clin_df <- as.data.frame(t_clin[-1, ], stringsAsFactors = FALSE)

# 提取关键信息：样本ID 和 分子亚型
# 注意：根据你的文件snippet，亚型所在的属性名是 "Integrated.Phenotype"
sample_metadata <- clin_df %>%
  dplyr:: select(sample_id = attrib_name, phenotype = Integrated.Phenotype) %>%
  filter(!is.na(phenotype) & phenotype != "" & phenotype != "NA") %>% # 去除无分类的样本
  filter(phenotype %in% c("Epithelial", "EMT", "Hypermutated"))       # 确保只保留这三类

# 打印一下检查
message("成功提取元数据，样本数：", nrow(sample_metadata))
table(sample_metadata$phenotype)

plot_integrated_heatmap <- function(expr_matrix, meta_df, platform_name, is_log_scale = FALSE) {
  
  # 1. 样本对齐
  # 找出表达矩阵和元数据共有的样本 ID
  common_samples <- intersect(colnames(expr_matrix), meta_df$sample_id)
  
  if(length(common_samples) < 3) {
    message("警告：", platform_name, " 匹配到的样本过少，跳过绘图。")
    return(NULL)
  }
  
  # 提取共有样本的数据
  dat <- expr_matrix[, common_samples]
  meta <- meta_df[match(common_samples, meta_df$sample_id), ]
  
  # 2. 数据过滤与填充
  # 过滤缺失值过多的蛋白 (超过 50% NA 则剔除)
  dat <- dat[rowSums(is.na(dat)) < ncol(dat) * 0.5, ]
  
  # 简单的 NA 填充 (使用行均值，仅为了计算 P 值和热图展示)
  dat_imputed <- t(apply(dat, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  
  # 如果不是 Log 尺度 (如 LF Counts)，进行 log2(x+1) 处理
  if(!is_log_scale) {
    dat_imputed <- log2(dat_imputed + 1)
  }
  
  # 3. 筛选 Top 100 差异蛋白 (ANOVA)
  # 计算每个蛋白在三种亚型间的差异 P 值
  pvalues <- apply(dat_imputed, 1, function(x) {
    tryCatch({
      # 排除极少数方差为0的情况
      if(var(x) == 0) return(1)
      res <- aov(x ~ as.factor(meta$phenotype))
      summary(res)[[1]][["Pr(>F)"]][1]
    }, error = function(e) 1)
  })
  
  # 排序并取前 100
  top100_ids <- names(sort(pvalues)[1:100])
  dat_top <- dat_imputed[top100_ids, ]
  
  # 4. Z-score 标准化 (Row scaling)
  dat_scaled <- t(scale(t(dat_top)))
  
  # 5. 构建热图注释
  # 使用特定的配色方案
  subtype_colors <- c("Epithelial" = "#66C2A5", 
                      "EMT" = "#FC8D62", 
                      "Hypermutated" = "#8DA0CB")
  
  ha <- HeatmapAnnotation(
    Subtype = meta$phenotype,
    col = list(Subtype = subtype_colors),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.4, "cm")
  )
  
  # 6. 绘制热图
  ht <- Heatmap(dat_scaled,
                name = "Expression\n(Z-score)",
                top_annotation = ha,
                col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")), # 蓝白红配色
                show_column_names = FALSE,  # 样本名不显示 (太挤)
                show_row_names = FALSE,     # 蛋白名不显示 (如果要看具体蛋白改为 TRUE)
                cluster_columns = TRUE,     # 对列聚类
                cluster_rows = TRUE,        # 对行聚类
                column_split = meta$phenotype, # 强制按亚型拆分
                column_title = paste0(platform_name, " Top 100 DE Proteins"),
                row_dend_width = unit(2, "cm")
  )
  
  return(ht)
}
# --- 绘制 TMT 平台 (Log Ratio, 不需要额外 Log) ---
pdf("Figure_Heatmap_Top100_TMT.pdf", width = 10, height = 8)
ht_tmt <- plot_integrated_heatmap(tmt_raw, sample_metadata, "TMT Platform", is_log_scale = TRUE)
draw(ht_tmt, merge_legend = TRUE)
dev.off()

# --- 绘制 LF 平台 (Counts, 需要 Log) ---
pdf("Figure_Heatmap_Top100_LF.pdf", width = 10, height = 8)
ht_lf <- plot_integrated_heatmap(lf_raw, sample_metadata, "LF Platform", is_log_scale = FALSE)
draw(ht_lf, merge_legend = TRUE)
dev.off()

message("热图绘制完成！请查看 PDF 文件。")
##提取TOP100蛋白保存###
library(openxlsx)
library(dplyr)
library(tidyr)
get_top100_table <- function(expr_matrix, meta_df, platform_name, is_log_scale = FALSE) {
  
  # 1. 样本与元数据对齐
  common_samples <- intersect(colnames(expr_matrix), meta_df$sample_id)
  dat <- expr_matrix[, common_samples]
  meta <- meta_df[match(common_samples, meta_df$sample_id), ]
  
  # 2. 预处理
  dat <- dat[rowSums(is.na(dat)) < ncol(dat) * 0.5, ]
  dat_imputed <- t(apply(dat, 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }))
  if(!is_log_scale) dat_imputed <- log2(dat_imputed + 1)
  
  # 3. 计算 ANOVA P 值
  pvalues <- apply(dat_imputed, 1, function(x) {
    tryCatch({
      if(var(x) == 0) return(1)
      res <- aov(x ~ as.factor(meta$phenotype))
      summary(res)[[1]][["Pr(>F)"]][1]
    }, error = function(e) 1)
  })
  
  # 4. 构建详细的信息表
  # 计算各组的平均表达量 (Mean Expression per Subtype)
  group_means <- as.data.frame(t(apply(dat_imputed, 1, function(x) {
    tapply(x, meta$phenotype, mean, na.rm = TRUE)
  })))
  
  # 整合表格
  full_table <- data.frame(
    Protein_Symbol = names(pvalues),
    ANOVA_P_value = pvalues,
    group_means,
    stringsAsFactors = FALSE
  ) %>%
    # 修正 P 值，添加显著性标记
    mutate(FDR = p.adjust(ANOVA_P_value, method = "BH")) %>%
    # 按 P 值从低到高排序，取前 100
    arrange(ANOVA_P_value) %>%
    head(100)
  
  # 为表格添加平台来源标记
  full_table$Platform = platform_name
  
  return(full_table)
}
# --- 1. 获取两个平台的表格 ---
# 假设 tmt_raw, lf_raw, sample_metadata 已经在您的环境中
tmt_top100_df <- get_top100_table(tmt_raw, sample_metadata, "TMT", is_log_scale = TRUE)
lf_top100_df <- get_top100_table(lf_raw, sample_metadata, "LF", is_log_scale = FALSE)

# --- 2. 创建并保存 Excel ---
wb <- createWorkbook()

# 添加 TMT 工作表
addWorksheet(wb, "TMT_Top100_Proteins")
writeData(wb, "TMT_Top100_Proteins", tmt_top100_df)

# 添加 LF 工作表
addWorksheet(wb, "LF_Top100_Proteins")
writeData(wb, "LF_Top100_Proteins", lf_top100_df)

# 保存文件
saveWorkbook(wb, "Supplementary_Table_Top100_Proteins.xlsx", overwrite = TRUE)

message("表格已成功保存为：Supplementary_Table_Top100_Proteins.xlsx")
common_top_proteins <- intersect(tmt_top100_df$Protein_Symbol, lf_top100_df$Protein_Symbol)
print(common_top_proteins)