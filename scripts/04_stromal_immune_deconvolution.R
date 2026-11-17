install.packages("tidyestimate")
library(DESeq2)
library(tidyverse)
library(msigdbr)
library(tidyestimate)

# 1. Load clean counts and clinical metadata
counts <- readRDS("data/processed/clean_counts.rds")
clinical <- readRDS("data/processed/clean_clinical.rds")

# 2. Compute VST-normalized matrix for scoring
keep <- rowSums(counts >= 10) >= (ncol(counts) * 0.2)
dds_tmp <- DESeqDataSetFromMatrix(countData = counts[keep, ], colData = clinical, design = ~ 1)
vsd_tmp <- vst(dds_tmp, blind = TRUE)
vst_df <- as.data.frame(assay(vsd_tmp)) %>%
  rownames_to_column(var = "gene_id") %>%
  mutate(ensembl_clean = gsub("\\..*", "", gene_id))

# 3. Map Ensembl IDs to unique HGNC Gene Symbols
gene_dict <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(ensembl_gene, gene_symbol) %>%
  distinct()

vst_hgnc <- vst_df %>%
  inner_join(gene_dict, by = c("ensembl_clean" = "ensembl_gene")) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "") %>%
  dplyr::select(-gene_id, -ensembl_clean) %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), mean), .groups = "drop")

# 4. Run ESTIMATE scoring via tidyestimate (genes as rownames, tidy = FALSE)
message("Calculating Stromal, Immune, and ESTIMATE scores...")
vst_mat_hgnc <- vst_hgnc %>%
  column_to_rownames(var = "gene_symbol") %>%
  as.data.frame()

est_scores <- vst_mat_hgnc %>%
  filter_common_genes(id = "hgnc_symbol", tidy = FALSE, tell_missing = FALSE, find_alias = FALSE) %>%
  estimate_score(is_affymetrix = FALSE)

# 5. Merge scores with Subtype labels
scores_merged <- est_scores %>%
  left_join(clinical %>% dplyr::select(barcode, Subtype), by = c("sample" = "barcode")) %>%
  filter(!is.na(Subtype))

write.csv(scores_merged, "results/tables/estimate_microenvironment_scores.csv", row.names = FALSE)

# 6. Statistical comparison (Wilcoxon rank-sum test)
p_stromal <- wilcox.test(stromal ~ Subtype, data = scores_merged)$p.value
p_immune  <- wilcox.test(immune ~ Subtype, data = scores_merged)$p.value
p_est     <- wilcox.test(estimate ~ Subtype, data = scores_merged)$p.value

message(sprintf("Wilcoxon p-values -> Stromal: %.4g | Immune: %.4g | ESTIMATE: %.4g",
                p_stromal, p_immune, p_est))

# 7. Generate Figure 5: Microenvironment Boxplots
scores_long <- scores_merged %>%
  pivot_longer(cols = c(stromal, immune, estimate), names_to = "Score_Type", values_to = "Score") %>%
  mutate(Score_Type = factor(Score_Type,
                             levels = c("stromal", "immune", "estimate"),
                             labels = c("Stromal Score", "Immune Score", "Composite ESTIMATE")))

p_tme <- ggplot(scores_long, aes(x = Subtype, y = Score, fill = Subtype)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  facet_wrap(~ Score_Type, scales = "free_y") +
  scale_fill_manual(values = c("Classical" = "#2B5C8F", "Basal_like" = "#D95F02")) +
  labs(
    title = "Tumor Microenvironment Infiltration: Basal-like vs. Classical PDAC",
    subtitle = sprintf("Wilcoxon p-values: Stromal = %.2e | Immune = %.2e", p_stromal, p_immune),
    x = "Molecular Subtype",
    y = "ssGSEA Enrichment Score"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

ggsave("results/figures/figure5_microenvironment_scores.png", p_tme, width = 10, height = 5.5, dpi = 300)
message("Figure 5 saved to results/figures/figure5_microenvironment_scores.png")