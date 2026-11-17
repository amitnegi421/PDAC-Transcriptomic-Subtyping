library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

#1. query GDC for TCGA-PAAD raw STAR counts
message("Querying GDC API for TCGA-PAAD...")
query_paad <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

#2. Download and assemble into a SummarizedExperiment
GDCdownload(query_paad, method = "api", files.per.chunk = 10)
paad_se <- GDCprepare(query_paad)

#3. Extract metadata and raw count matrix
raw_counts <- assay(paad_se, "unstranded")
clinical_meta <- as.data.frame(colData(paad_se))

#4. Filter for primary solid tumors (code '01')
tumor_sample <- clinical_meta %>% 
  filter(sample_type == "Primary Tumor" | grepl("-01A", barcode))

#5.Biological Curation: Retain verified Pancreatic Ductal Adenocarcinoma(PDAC)
# Excluding neuroendocrine tumors and non-ductal variants

clean_clinical <- tumor_sample %>%
  filter(!is.na(days_to_death) | !is.na(days_to_last_follow_up)) %>%
  filter(!grepl("neuroendocrine|acinar|colloid", primary_diagnosis, ignore.case = TRUE))

#Synchronise count matrix with curated patients
filtered_counts <- raw_counts[, clean_clinical$barcode]
message(paste("Retained", ncol(filtered_counts), "curated PDAC tumor samples."))
 
#save clean working objects
saveRDS(filtered_counts, file = "data/processed/clean_counts.rds")
saveRDS(clean_clinical, file = "data/processed/clean_clinical.rds")
           