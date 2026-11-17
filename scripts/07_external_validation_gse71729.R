if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
library(msigdbr)
library(pheatmap)

# 1. Load your stromal-adjusted DESeq2 results from Script 06
res_adj <- read.csv("results/tables/deseq2_stromal_adjusted.csv")

gene_dict <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(ensembl_gene, gene_symbol) %>%
  distinct()

res_mapped <- res_adj %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(ensembl_clean = gsub("\\..*", "", gene_id)) %>%
  inner_join(gene_dict, by = c("ensembl_clean" = "ensembl_gene")) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")

# Top 25 Tumor-Intrinsic Basal-like and Classical markers
top_basal <- res_mapped %>%
  filter(padj < 0.01 & log2FoldChange > 1.5) %>%
  arrange(padj) %>% distinct(gene_symbol) %>% slice_head(n = 25) %>% pull(gene_symbol)

top_classical <- res_mapped %>%
  filter(padj < 0.01 & log2FoldChange < -1.5) %>%
  arrange(padj) %>% distinct(gene_symbol) %>% slice_head(n = 25) %>% pull(gene_symbol)

# 2. Download independent validation cohort GSE71729 from NCBI GEO
message("Downloading GSE71729 from NCBI GEO...")
gse <- getGEO(filename = "data/raw/GSE71729_series_matrix.txt", getGPL = FALSE)

expr_mat <- exprs(gse)
pdata <- pData(gse)
fdata <- fData(gse)
head(rownames(expr_mat))
# Inspect metadata columns in console
message("GSE71729 metadata columns:")
print(colnames(pdata))

# 3. Ensure expression matrix row names are HGNC gene symbols
if (!any(c("GATA6", "KRT17") %in% rownames(expr_mat))) {
  sym_col <- grep("symbol", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sym_col)) {
    rownames(expr_mat) <- make.unique(as.character(fdata[[sym_col]]))
  }
}

# 4. Score GSE71729 samples using our TCGA-derived Tumor-Intrinsic Signature
basal_present <- intersect(top_basal, rownames(expr_mat))
class_present <- intersect(top_classical, rownames(expr_mat))

message(sprintf("Matched %d/25 Basal and %d/25 Classical intrinsic genes in GSE71729.",
                length(basal_present), length(class_present)))

# Z-score genes across GSE71729 samples and compute Intrinsic Subtype Score
expr_sub <- expr_mat[c(basal_present, class_present), , drop = FALSE]
expr_z <- t(scale(t(expr_sub)))

basal_score <- colMeans(expr_z[basal_present, , drop = FALSE], na.rm = TRUE)
class_score <- colMeans(expr_z[class_present, , drop = FALSE], na.rm = TRUE)
intrinsic_score <- basal_score - class_score

val_df <- data.frame(
  Sample = colnames(expr_mat),
  Basal_Score = basal_score,
  Classical_Score = class_score,
  Intrinsic_Score = intrinsic_score,
  Predicted_Subtype = ifelse(intrinsic_score > 0, "Basal_like", "Classical")
)

write.csv(val_df, "results/tables/gse71729_validation_scores.csv", row.names = FALSE)
print(table(val_df$Predicted_Subtype))

# 5. Generate Figure 8: Validation Heatmap on GSE71729
col_order <- order(val_df$Intrinsic_Score)
ann_col <- data.frame(Predicted_Subtype = val_df$Predicted_Subtype)
rownames(ann_col) <- val_df$Sample

pheatmap(
  expr_z[, col_order],
  annotation_col = ann_col,
  annotation_colors = list(Predicted_Subtype = c("Classical" = "#2B5C8F", "Basal_like" = "#D95F02")),
  cluster_cols = FALSE,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "External Validation (GSE71729): Tumor-Intrinsic Basal vs. Classical Signature",
  filename = "results/figures/figure8_external_validation_gse71729.png",
  width = 10,
  height = 8
)
message("Figure 8 saved to results/figures/figure8_external_validation_gse71729.png")

