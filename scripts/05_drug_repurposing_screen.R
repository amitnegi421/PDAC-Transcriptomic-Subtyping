install.packages("enrichR")
library(tidyverse)
library(msigdbr)
library(enrichR)

# 1. Load DESeq2 results and map Ensembl IDs to HGNC symbols
res_df <- read.csv("results/tables/deseq2_basal_vs_classical.csv")
gene_dict <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(ensembl_gene, gene_symbol) %>%
  distinct()

res_mapped <- res_df %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(ensembl_clean = gsub("\\..*", "", gene_id)) %>%
  inner_join(gene_dict, by = c("ensembl_clean" = "ensembl_gene")) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")

# 2. Extract top 200 Basal-upregulated and top 200 Basal-downregulated genes
basal_up <- res_mapped %>%
  filter(padj < 0.01 & log2FoldChange > 1.5) %>%
  arrange(padj) %>%
  distinct(gene_symbol) %>%
  slice_head(n = 200) %>%
  pull(gene_symbol)

basal_down <- res_mapped %>%
  filter(padj < 0.01 & log2FoldChange < -1.5) %>%
  arrange(padj) %>%
  distinct(gene_symbol) %>%
  slice_head(n = 200) %>%
  pull(gene_symbol)

message(sprintf("Querying LINCS L1000 with %d Basal-up and %d Basal-down genes...",
                length(basal_up), length(basal_down)))

# 3. Query Enrichr LINCS L1000 Chemical Perturbation libraries
setEnrichrSite("Enrichr")
res_rev_up   <- enrichr(basal_up,   "LINCS_L1000_Chem_Pert_down")[["LINCS_L1000_Chem_Pert_down"]]
res_rev_down <- enrichr(basal_down, "LINCS_L1000_Chem_Pert_up")[["LINCS_L1000_Chem_Pert_up"]]

# 4. Tag direction and combine top candidates
top_candidates <- bind_rows(
  res_rev_up   %>% mutate(Direction = "Downregulates Basal-Up Genes"),
  res_rev_down %>% mutate(Direction = "Upregulates Classical Genes")
) %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(desc(Combined.Score))

write.csv(top_candidates, "results/tables/lincs_drug_candidates.csv", row.names = FALSE)
message(paste("Saved", nrow(top_candidates), "significant drug-perturbation hits."))

# 5. Generate Figure 6: Top 15 Reversal Candidates
plot_df <- top_candidates %>%
  slice_head(n = 15) %>%
  mutate(Term_clean = str_trunc(Term, 45))

p_drugs <- ggplot(plot_df, aes(x = reorder(Term_clean, Combined.Score), y = Combined.Score, fill = Direction)) +
  geom_col(alpha = 0.9) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Downregulates Basal-Up Genes" = "#D95F02",
    "Upregulates Classical Genes"  = "#2B5C8F"
  )) +
  labs(
    title = "Top In Silico Drug Repurposing Candidates (LINCS L1000)",
    subtitle = "Compounds reversing the Basal-like PDAC transcriptional signature (FDR < 0.05)",
    x = "Compound / Perturbation",
    y = "Enrichr Combined Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("results/figures/figure6_drug_repurposing.png", p_drugs, width = 10, height = 6.5, dpi = 300)
message("Figure 6 saved to results/figures/figure6_drug_repurposing.png")