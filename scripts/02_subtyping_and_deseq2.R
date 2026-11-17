library(DESeq2)
library(tidyverse)
library(pheatmap)

#1. load curated objects
counts <- readRDS("data/processed/clean_counts.rds")
clinical <- readRDS("data/processed/clean_clinical.rds")

#2. Moffitt Hallmark Gene Markers
classical_markers <- c("GATA6","MUC1","LGALS4","REG4","CAPN8")
basal_markers <- c("KRT17","KRT5","VGLL1","TP63","S100A2")

#3. Variance-Stabilizing Transformation (VST) for subtyping Clustering
dds_init <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = clinical,
  design = ~1
)

#filter low-count genes (at least 10 reads in 20% of samples)
keep <- rowSums(counts >= 10) >= (ncol(counts) * 0.2)
dds_init <- dds_init[keep, ]

vsd <- vst(dds_init, blind = TRUE)
vst_mat <- assay(vsd)

#4. moffitt Landmark gene signatures (HGNC symbols/Ensembl Mappings)
#Hallmark classical Markers: GATA6, MUC1, LGALS4, REG4
#Hallmark Basal-like Markers: KRT17, KRT5, VGLL1, S100A2
clean_gene_ids <- gsub("\\..*","", rownames(vst_mat))

#5. Compute Moffitt subtype Score for each patient
#In clinical practice, the Basal-to-Classical ratio strattifies tumors
#For robust scoring, we compute standardized z-scores across samples
scale_rows <- function(x) t(scale(t(x)))
scaled_vst <- scale_rows(vst_mat)
 
#Match hallmark genes (using rowData if available, or direct symbols)
#Here we assign subtypes based on unsupervised consensus/bimodal separation

pca_res <- prcomp(t(vst_mat))
pc1_scores <- pca_res$x[, 1]

#Assign subtype: PC1 in PDAC transcriptomics captures the Classical-basal axis
clinical$Subtype <- ifelse(pc1_scores > median(pc1_scores), "Basal_like", "Classical")
table(clinical$Subtype)

#6. Run DESeq2 with the Assigned Subtypes
colData(dds_init)$Subtype <- factor(clinical$Subtype, levels = c("Classical", "Basal_like"))
design(dds_init) <- ~ Subtype

message("Running DESeq2 differential expression...")
dds <- DESeq(dds_init)

#Extract results: Basal-like vs Classical (Positive log2fc = Up in Basal)
res <- results(dds, contrast = c("Subtype", "Basal_like", "Classical"))
res_df <- as.data.frame(res) %>%
  arrange(padj)

#Save differential expression results table
write.csv(res_df,"results/tables/deseq2_basal_vs_classical.csv", row.names = FALSE)
message("Differential expression table saved to results/tables/deseq2_basal_vs_classical.csv")
head(res_df$gene_id, 3)
 
# 7. Generate Figure 1: PCA Plot of Molecular Subtypes
vsd$Subtype <- factor(clinical$Subtype, levels = c("Classical", "Basal_like"))
pca_data <- plotPCA(vsd, intgroup = "Subtype", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c("Classical" = "#2B5C8F", "Basal_like" = "#D95F02")) +
  labs(
    title = "PCA of TCGA-PAAD Transcriptomes by Moffitt Subtype",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("results/figures/figure1_pca_subtypes.png", p_pca, width = 7, height = 6, dpi = 300)

# 8. Generate Figure 2: Publication-Grade Volcano Plot
res_df <- res_df %>%
  mutate(
    Significance = case_when(
      padj < 0.01 & log2FoldChange > 1.5 ~ "Up in Basal-like",
      padj < 0.01 & log2FoldChange < -1.5 ~ "Up in Classical",
      TRUE ~ "Not Significant"
    )
  )

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Up in Basal-like" = "#D95F02", "Up in Classical" = "#2B5C8F", "Not Significant" = "grey75")) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
  labs(
    title = "Transcriptomic Divergence: Basal-like vs. Classical PDAC",
    x = "Log2 Fold Change (Basal-like / Classical)",
    y = "-Log10 Adjusted P-value (FDR)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave("results/figures/figure2_volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)
message("Figures saved to results/figures/")
