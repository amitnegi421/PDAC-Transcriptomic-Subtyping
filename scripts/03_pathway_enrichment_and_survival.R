if (!requireNamespace("msigdbr", quietly = TRUE)) 
  install.packages("msigdbr")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
library(DESeq2)
library(tidyverse)
library(survival)
library(survminer)
library(clusterProfiler)
library(msigdbr)

# 1. Load curated clinical data
clinical <- readRDS("data/processed/clean_clinical.rds")
if (!"Subtype" %in% colnames(clinical)) {
  message("Subtype column missing — computing from clean_counts.rds and saving...")
  counts <- readRDS("data/processed/clean_counts.rds")
  keep <- rowSums(counts >= 10) >= (ncol(counts) * 0.2)
  dds_tmp <- DESeqDataSetFromMatrix(countData = counts[keep, ], colData = clinical, design = ~ 1)
  vsd_tmp <- vst(dds_tmp, blind = TRUE)
  pc1_scores <- prcomp(t(assay(vsd_tmp)))$x[, 1]
  clinical$Subtype <- factor(
    ifelse(pc1_scores > median(pc1_scores), "Basal_like", "Classical"),
    levels = c("Classical", "Basal_like")
  )
  saveRDS(clinical, "data/processed/clean_clinical.rds")
  message("Updated clean_clinical.rds with Subtype.")
}
# Ensure deseq2_basal_vs_classical.csv has gene_id as column 1
res_path <- "results/tables/deseq2_basal_vs_classical.csv"
res_df <- if (file.exists(res_path)) read.csv(res_path) else NULL

if (is.null(res_df) || !"gene_id" %in% colnames(res_df)) {
  message("Extracting DESeq2 results with gene_id column...")
  if (!exists("dds")) {
    counts <- readRDS("data/processed/clean_counts.rds")
    keep <- rowSums(counts >= 10) >= (ncol(counts) * 0.2)
    dds <- DESeqDataSetFromMatrix(countData = counts[keep, ], colData = clinical, design = ~ Subtype)
    dds <- DESeq(dds)
  }
  res <- results(dds, contrast = c("Subtype", "Basal_like", "Classical"))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  write.csv(res_df, res_path, row.names = FALSE)
}

# ==============================================================================
# PART A: Kaplan-Meier Overall Survival Modeling (Figure 3)
# ==============================================================================
clinical_surv <- clinical %>%
  mutate(
    time_days = ifelse(!is.na(days_to_death), as.numeric(days_to_death), as.numeric(days_to_last_follow_up)),
    time_months = time_days / 30.4375,
    status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(time_months) & time_months > 0)

fit_subtype <- survfit(Surv(time_months, status) ~ Subtype, data = clinical_surv)

km_plot <- ggsurvplot(
  fit_subtype,
  data = clinical_surv,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("#D95F02", "#2B5C8F"),
  legend.title = "Molecular Subtype",
  legend.labs = c("Basal-like", "Classical"),
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival in TCGA-PAAD: Basal-like vs. Classical",
  ggtheme = theme_classic(base_size = 14)
)

png("results/figures/figure3_kaplan_meier_survival.png", width = 8, height = 7, units = "in", res = 300)
print(km_plot)
dev.off()
message("Figure 3 saved to results/figures/figure3_kaplan_meier_survival.png")

# ==============================================================================
# PART B: Hallmark Gene Set Enrichment Analysis (Figure 4)
# ==============================================================================
gene_dict <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(ensembl_gene, gene_symbol) %>%
  distinct()

res_clean <- res_df %>%
  filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(pvalue)) %>%
  mutate(
    ensembl_clean = gsub("\\..*", "", gene_id),
    rank_metric = sign(log2FoldChange) * -log10(pvalue + 1e-300)
  ) %>%
  inner_join(gene_dict, by = c("ensembl_clean" = "ensembl_gene")) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")

res_unique <- res_clean %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_metric))

message(paste("Mapped", nrow(res_unique), "unique gene symbols for GSEA."))

gene_ranks <- setNames(res_unique$rank_metric, res_unique$gene_symbol)
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_results <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_df <- as.data.frame(gsea_results)
write.csv(gsea_df, "results/tables/gsea_hallmark_results.csv", row.names = FALSE)
message(paste("Found", nrow(gsea_df), "enriched Hallmark pathways."))

if (nrow(gsea_df) > 0) {
  p_gsea <- dotplot(gsea_results, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign) +
    scale_color_gradient(low = "#D95F02", high = "#2B5C8F") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Hallmark Pathways Enriched in Basal-like vs. Classical PDAC",
      x = "Gene Ratio",
      y = "Pathway"
    ) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("results/figures/figure4_gsea_hallmarks.png", p_gsea, width = 10, height = 7, dpi = 300)
  message("Figure 4 saved to results/figures/figure4_gsea_hallmarks.png")
}