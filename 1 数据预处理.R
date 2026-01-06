library(tidyverse)
library(limma)
library(impute) 
setwd("C:/Users/zcy15/Downloads/Phenome homework/data")
# ==============================================================================
# 1. 读取数据
# ==============================================================================
# 读取临床数据
clinical <- read_tsv("Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi") %>%
  # 假设第一列是 attrib_name，如果文件结构不同需调整
  column_to_rownames(var = colnames(.)[1]) %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  select(Sample, Integrated.Phenotype) %>%
  filter(!is.na(Integrated.Phenotype), 
         Integrated.Phenotype %in% c("Hypermutated", "Epithelial", "EMT"))

# 读取 Label-Free (LF)
lf_raw <- read_tsv("Human__CPTAC_COAD__VU__Proteome__QExact__03_01_2017__BCM__Gene__VU_Tumor_LF_UnsharedCounts.cct") %>%
  column_to_rownames(var = colnames(.)[1])

# 读取 TMT
tmt_raw <- read_tsv("Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct") %>%
  column_to_rownames(var = colnames(.)[1])

# ==============================================================================
# 2. 修正后的预处理函数
# ==============================================================================
preprocess_proteome <- function(prot_df, clinical_df, type = "LF", min_valid_per_group = 0.5) {
  
  # [步骤A] 数据清洗与对齐
  # 1. 找到公共样本
  common_samples <- intersect(colnames(prot_df), clinical_df$Sample)
  message("匹配到的样本数: ", length(common_samples))
  
  # 2. 对齐矩阵
  prot_df <- prot_df[, common_samples]
  clin_matched <- clinical_df %>% filter(Sample %in% common_samples) %>% arrange(match(Sample, common_samples))
  
  # 3. 提取表型因子
  pheno <- factor(clin_matched$Integrated.Phenotype, levels = c("Epithelial", "EMT", "Hypermutated"))
  
  # [步骤B] 缺失值处理
  expr_mat <- as.matrix(prot_df)
  
  # *修正1*: Label-free 的 0 需要转为 NA
  if(type == "LF") {
    expr_mat[expr_mat == 0] <- NA
  }
  
  # *修正2*: 过滤逻辑 (按组计算)
  # 创建一个逻辑矩阵，判断每个值是否有效
  is_valid <- !is.na(expr_mat)
  keep_genes <- rep(FALSE, nrow(expr_mat))
  
  # 循环每个组，只要有一个组满足条件就保留 (Union)
  for(g in levels(pheno)) {
    samples_in_group <- which(pheno == g)
    # 计算该组内每个蛋白的有效比例
    valid_rate <- rowMeans(is_valid[, samples_in_group, drop=FALSE])
    keep_genes <- keep_genes | (valid_rate >= min_valid_per_group)
  }
  
  expr_filtered <- expr_mat[keep_genes, ]
  message("过滤后保留蛋白数: ", nrow(expr_filtered))
  
  # [步骤C] 缺失值填补
  if(nrow(expr_filtered) > 0) {
    if(type == "LF") {
      # Label-free: 使用最小值/2 (模拟低丰度)
      message("使用 Min/2 填补 (Label-Free)")
      min_val <- min(expr_filtered, na.rm = TRUE) # 全局最小或行最小均可，这里用每行最小更精细
      
      # 按行填补
      expr_imputed <- t(apply(expr_filtered, 1, function(x) {
        x[is.na(x)] <- min(x, na.rm=TRUE) / 2
        return(x)
      }))
      
    } else {
      # TMT: 推荐使用 KNN (因为是 LogRatio，直接补最小值会产生偏差)
      message("使用 KNN 填补 (TMT)")
      # impute.knn 需要 (行=基因, 列=样本)
      imputed_res <- impute.knn(expr_filtered, k = 10)
      expr_imputed <- imputed_res$data
    }
  } else {
    expr_imputed <- expr_filtered
  }
  
  # [步骤D] 数据转换 (LF 需要 Log2 转换用于 limma)
  if(type == "LF") {
    expr_imputed <- log2(expr_imputed + 1) # +1 防止 log(0)
  }
  
  # *修正3*: 返回列表，确保 expr 是 (Protein x Sample)
  list(expr = expr_imputed, pheno = pheno, samples = common_samples)
}

# ==============================================================================
# 3. 执行预处理
# ==============================================================================
# 注意：LF 数据我加了 log2 转换步骤，因为 limma 也是基于正态分布假设
lf_processed  <- preprocess_proteome(lf_raw, clinical, type = "LF", min_valid_per_group = 0.5)
tmt_processed <- preprocess_proteome(tmt_raw, clinical, type = "TMT", min_valid_per_group = 0.5)
# 推荐方案2 - 样本×基因格式
expr_transposed <- t(lf_processed$expr)
colnames(expr_transposed) <- rownames(lf_processed$expr)
rownames(expr_transposed) <- lf_processed$samples

# 创建包含表型的完整表
full_data <- data.frame(
  sample_id = lf_processed$samples,
  phenotype = lf_processed$pheno,
  expr_transposed,
  check.names = FALSE
)

# 保存
write.csv(full_data, "LF_analysis_ready.csv", row.names = FALSE)
#####
expr_transposed <- t(tmt_processed$expr)
colnames(expr_transposed) <- rownames(tmt_processed$expr)
rownames(expr_transposed) <- tmt_processed$samples

# 创建包含表型的完整表
full_data <- data.frame(
  sample_id = tmt_processed$samples,
  phenotype = tmt_processed$pheno,
  expr_transposed,
  check.names = FALSE
)

# 保存
write.csv(full_data, "tmt_analysis_ready.csv", row.names = FALSE)

# ==============================================================================
# 4. 差异表达分析函数 (Limma) - 修正 coef 提取命名
# ==============================================================================
run_diff_analysis <- function(processed_data) {
  expr <- processed_data$expr
  pheno <- processed_data$pheno
  
  # 设计矩阵
  design <- model.matrix(~ 0 + pheno)
  colnames(design) <- levels(pheno)
  
  fit <- lmFit(expr, design)
  
  # 定义对比，使用你希望的名称
  contrast_matrix <- makeContrasts(
    EMT_vs_Epithelial   = EMT - Epithelial,
    Hyper_vs_Epithelial = Hypermutated - Epithelial,
    Hyper_vs_EMT        = Hypermutated - EMT,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  print(colnames(fit2$coefficients))  
  # 提取结果时，使用 makeContrasts 中定义的名称
  results <- list(
    EMT_vs_Epithelial = topTable(fit2, coef = "EMT_vs_Epithelial",   number = Inf) %>% rownames_to_column("Protein"),
    Hyper_vs_Epithelial = topTable(fit2, coef = "Hyper_vs_Epithelial", number = Inf) %>% rownames_to_column("Protein"),
    Hyper_vs_EMT = topTable(fit2, coef = "Hyper_vs_EMT", number = Inf) %>% rownames_to_column("Protein")
  )
  
  return(results)
}
# 在 fit2 <- eBayes(fit2) 之后添加
print(colnames(fit2$coefficients))
# ==============================================================================
# 5. 运行分析并保存 - 统一使用正确的 key
# ==============================================================================
lf_diff <- run_diff_analysis(lf_processed)
tmt_diff <- run_diff_analysis(tmt_processed)

# 查看显著蛋白数量（可选）
print("LF - EMT vs Epithelial:")
print(table(lf_diff$EMT_vs_Epithelial$adj.P.Val < 0.05))
print("TMT - EMT vs Epithelial:")
print(table(tmt_diff$EMT_vs_Epithelial$adj.P.Val < 0.05))

print("LF - Hyper vs EMT:")
print(table(lf_diff$Hyper_vs_EMT$adj.P.Val < 0.05))
print("TMT - Hyper vs EMT:")
print(table(tmt_diff$Hyper_vs_EMT$adj.P.Val < 0.05))

print("LF - Hyper vs Epithelial:")
print(table(lf_diff$Hyper_vs_Epithelial$adj.P.Val < 0.05))
print("TMT - Hyper vs Epithelial:")
print(table(tmt_diff$Hyper_vs_Epithelial$adj.P.Val < 0.05))

# 保存文件
write_csv(lf_diff$EMT_vs_Epithelial,   "LF_EMT_vs_Epithelial_Corrected.csv")
write_csv(tmt_diff$EMT_vs_Epithelial,  "tmt_EMT_vs_Epithelial_Corrected.csv")

write_csv(lf_diff$Hyper_vs_EMT,        "LF_Hyper_vs_EMT_Corrected.csv")
write_csv(tmt_diff$Hyper_vs_EMT,       "tmt_Hyper_vs_EMT_Corrected.csv")

write_csv(lf_diff$Hyper_vs_Epithelial, "LF_Hyper_vs_Epithelial_Corrected.csv")
write_csv(tmt_diff$Hyper_vs_Epithelial,"tmt_Hyper_vs_Epithelial_Corrected.csv")