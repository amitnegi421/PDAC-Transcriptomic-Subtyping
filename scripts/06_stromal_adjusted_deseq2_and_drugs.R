library(DESeq2)
library(tidyverse)
library(msigdbr)
library(enrichR)

# 1. Load clean counts, clinical metadata, and ESTIMATE scores
counts <- readRDS("data/processed/clean_counts.rds")
clinical <- readRDS("data/processed/clean_clinical.rds")
est_scores <- read.csv("results/tables/estimate_microenvironment_scores.csv")
res_unadj <- read.csv("results/tables/deseq2_basal_vs_classical.csv")

# 2. Merge stromal scores and standardize (z-score) for DESeq2 GLM stability
clinical_adj <- clinical %>%
  left_join(est_scores %>% dplyr::select(sample, stromal), by = c("barcode" = "sample")) %>%
  filter(!is.na(stromal)) %>%
  mutate(
    stromal_z = as.numeric(scale(stromal)),
    Subtype = factor(Subtype, levels = c("Classical", "Basal_like"))
  )

counts_adj <- counts[, clinical_adj$barcode]
keep <- rowSums(counts_adj >= 10) >= (ncol(counts_adj) * 0.2)

# 3. Fit Stromal-Adjusted DESeq2 Model (~ stromal_z + Subtype)
message("Fitting stromal-adjusted DESeq2 model (~ stromal_z + Subtype)...")
dds_adj <- DESeqDataSetFromMatrix(
  countData = counts_adj[keep, ],
  colData = clinical_adj,
  design = ~ stromal_z + Subtype
)
dds_adj <- DESeq(dds_adj)

res_adj <- results(dds_adj, contrast = c("Subtype", "Basal_like", "Classical"))
res_adj_df <- as.data.frame(res_adj)
res_adj_df$gene_id <- rownames(res_adj_df)
res_adj_df <- res_adj_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.csv(res_adj_df, "results/tables/deseq2_stromal_adjusted.csv", row.names = FALSE)

# 4. Quantify Stromal Confounding vs. Tumor-Intrinsic Drivers
comp_df <- res_unadj %>%
  dplyr::select(gene_id, lfc_unadj = log2FoldChange, padj_unadj = padj) %>%
  inner_join(
    res_adj_df %>% dplyr::select(gene_id, lfc_adj = log2FoldChange, padj_adj = padj),
    by = "gene_id"
  ) %>%
  filter(!is.na(padj_unadj) & !is.na(padj_adj)) %>%
  mutate(
    Category = case_when(
      padj_unadj < 0.01 & abs(lfc_unadj) > 1.5 & padj_adj < 0.01 & abs(lfc_adj) > 1.5 ~ "Tumor-Intrinsic Driver",
      padj_unadj < 0.01 & abs(lfc_unadj) > 1.5 ~ "Stromal-Confounded Artifact",
      padj_adj < 0.01 & abs(lfc_adj) > 1.5 ~ "Unmasked Intrinsic Gene",
      TRUE ~ "Not Significant"
    )
  )

print(table(comp_df$Category))

# 5. Generate Figure 7: Unadjusted vs. Stromal-Adjusted Fold Changes
p_comp <- ggplot(comp_df, aes(x = lfc_unadj, y = lfc_adj, color = Category)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c(
    "Tumor-Intrinsic Driver"      = "#2B5C8F",
    "Stromal-Confounded Artifact" = "#D95F02",
    "Unmasked Intrinsic Gene"     = "#1B9E77",
    "Not Significant"             = "grey82"
  )) +
  labs(
    title = "Isolating Tumor-Intrinsic Drivers via Stromal Adjustment",
    subtitle = "Comparison of Unadjusted (~ Subtype) vs. Adjusted (~ stromal_z + Subtype) Log2FC",
    x = "Unadjusted Log2 Fold Change",
    y = "Stromal-Adjusted Log2 Fold Change"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("results/figures/figure7_stromal_confounding_comparison.png", p_comp, width = 8.5, height = 6.5, dpi = 300)
message("Figure 7 saved to results/figures/figure7_stromal_confounding_comparison.png")

# 6. Re-screen LINCS L1000 using ONLY Tumor-Intrinsic Genes
gene_dict <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(ensembl_gene, gene_symbol) %>%
  distinct()

intrinsic_mapped <- res_adj_df %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(ensembl_clean = gsub("\\..*", "", gene_id)) %>%
  inner_join(gene_dict, by = c("ensembl_clean" = "ensembl_gene")) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")

intrinsic_up <- intrinsic_mapped %>%
  filter(padj < 0.01 & log2FoldChange > 1.5) %>%
  arrange(padj) %>% distinct(gene_symbol) %>% slice_head(n = 200) %>% pull(gene_symbol)

intrinsic_down <- intrinsic_mapped %>%
  filter(padj < 0.01 & log2FoldChange < -1.5) %>%
  arrange(padj) %>% distinct(gene_symbol) %>% slice_head(n = 200) %>% pull(gene_symbol)

setEnrichrSite("Enrichr")
drugs_up   <- enrichr(intrinsic_up,   "LINCS_L1000_Chem_Pert_down")[["LINCS_L1000_Chem_Pert_down"]]
drugs_down <- enrichr(intrinsic_down, "LINCS_L1000_Chem_Pert_up")[["LINCS_L1000_Chem_Pert_up"]]

intrinsic_drugs <- bind_rows(
  drugs_up   %>% mutate(Direction = "Downregulates Intrinsic Basal-Up"),
  drugs_down %>% mutate(Direction = "Upregulates Intrinsic Classical")
) %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(desc(Combined.Score))

write.csv(intrinsic_drugs, "results/tables/lincs_drug_candidates_intrinsic.csv", row.names = FALSE)
message(paste("Saved", nrow(intrinsic_drugs), "tumor-intrinsic drug candidates!"))