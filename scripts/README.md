# Integrative Transcriptomic Subtyping & In Silico Drug Repurposing in PDAC

A reproducible **Nextflow DSL2** and **R/Bioconductor** pipeline that stratifies Pancreatic Ductal Adenocarcinoma (PDAC) tumors into Classical and Basal-like subtypes, adjusts for stromal infiltration to isolate tumor-intrinsic drivers, validates the signature in an independent GEO cohort (`GSE71729`), and screens LINCS L1000 compounds that reverse the aggressive Basal-like phenotype.

## Pipeline Architecture
1. `01_fetch_and_curate_data.R` — Queries NCI GDC (`TCGAbiolinks`), filters for primary ductal adenocarcinomas, and exports raw STAR counts.
2. `02_subtyping_and_deseq2.R` — Applies VST normalization, assigns Classical vs. Basal-like subtypes, runs unadjusted `DESeq2`, and plots PCA & Volcano figures (Figures 1–2).
3. `03_pathway_enrichment_and_survival.R` — Fits Kaplan-Meier overall survival curves and runs MSigDB Hallmark `GSEA` (Figures 3–4).
4. `04_stromal_immune_deconvolution.R` — Quantifies stromal and immune infiltration via `tidyestimate` (Figure 5).
5. `05_drug_repurposing_screen.R` — Screens `LINCS_L1000_Chem_Pert` libraries via `enrichR` on unadjusted DEGs (Figure 6).
6. `06_stromal_adjusted_deseq2_and_drugs.R` — Fits a stromal-adjusted model (`~ stromal_z + Subtype`) to separate tumor-intrinsic drivers from stromal artifacts and re-screens LINCS L1000 (Figure 7).
7. `07_external_validation_gse71729.R` — Validates the 50-gene tumor-intrinsic signature in the independent `GSE71729` PDAC cohort (Figure 8).

## Quick Start
```bash
# Build container and run the full pipeline
docker build -t pdac-subtyping:latest .
nextflow run main.nf -with-docker pdac-subtyping:latest