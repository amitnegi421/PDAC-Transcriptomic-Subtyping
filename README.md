# PDAC-Transcriptomic-Subtyping
Reproducible Nextflow &amp; R/Bioconductor pipeline for stromal-adjusted PDAC subtyping and LINCS L1000 drug repurposing.
## Data Sources & Citations
- **Discovery Cohort (TCGA-PAAD):** [NCI Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/projects/TCGA-PAAD) — Cancer Genome Atlas Research Network (*Cancer Cell*, 2017).
- **Validation Cohort (GSE71729):** [NCBI GEO GSE71729](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729) — Moffitt et al., *Virtual microdissection identifies distinct tumor- and stroma-specific subtypes of pancreatic ductal adenocarcinoma* (*Nature Genetics*, 2015).
- **Tumor Microenvironment Scoring:** [ESTIMATE](https://www.nature.com/articles/ncomms3612) (Yoshihara et al., *Nature Communications*, 2013) implemented via [`tidyestimate`](https://cran.r-project.org/package=tidyestimate).
- **Drug Repurposing Screen:** [Enrichr](https://maayanlab.cloud/enrichr/) (`LINCS_L1000_Chem_Pert_down` and `LINCS_L1000_Chem_Pert_up` libraries).
