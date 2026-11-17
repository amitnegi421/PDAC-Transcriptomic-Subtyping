nextflow.enable.dsl=2

params.outdir = "results"

process FETCH_AND_CURATE {
    output:
    path "data/processed/clean_counts.rds", emit: counts
    path "data/processed/clean_clinical.rds", emit: clinical

    script:
    """
    mkdir -p data/processed
    Rscript ${projectDir}/scripts/01_fetch_and_curate_data.R
    """
}

process SUBTYPING_AND_DESEQ2 {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path clinical

    output:
    path "results/tables/deseq2_basal_vs_classical.csv", emit: deg_table
    path "data/processed/clean_clinical.rds", emit: clinical_subtyped
    path "results/figures/*.png"

    script:
    """
    mkdir -p data/processed results/tables results/figures
    Rscript ${projectDir}/scripts/02_subtyping_and_deseq2.R
    """
}

process PATHWAY_AND_SURVIVAL {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path clinical_subtyped
    path deg_table

    output:
    path "results/tables/gsea_hallmark_results.csv"
    path "results/figures/*.png"

    script:
    """
    mkdir -p results/tables results/figures
    Rscript ${projectDir}/scripts/03_pathway_enrichment_and_survival.R
    """
}

process TME_DECONVOLUTION {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path clinical_subtyped

    output:
    path "results/tables/estimate_microenvironment_scores.csv", emit: est_scores
    path "results/figures/*.png"

    script:
    """
    mkdir -p results/tables results/figures
    Rscript ${projectDir}/scripts/04_stromal_immune_deconvolution.R
    """
}

process DRUG_REPURPOSING {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path deg_table

    output:
    path "results/tables/lincs_drug_candidates.csv"
    path "results/figures/*.png"

    script:
    """
    mkdir -p results/tables results/figures
    Rscript ${projectDir}/scripts/05_drug_repurposing_screen.R
    """
}

process STROMAL_ADJUSTED_DESEQ2_AND_DRUGS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path clinical_subtyped
    path deg_table
    path est_scores

    output:
    path "results/tables/deseq2_stromal_adjusted.csv", emit: deg_adj
    path "results/tables/lincs_drug_candidates_intrinsic.csv"
    path "results/figures/*.png"

    script:
    """
    mkdir -p results/tables results/figures
    Rscript ${projectDir}/scripts/06_stromal_adjusted_deseq2_and_drugs.R
    """
}

process EXTERNAL_VALIDATION {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path deg_adj

    output:
    path "results/tables/gse71729_validation_scores.csv"
    path "results/figures/*.png"

    script:
    """
    mkdir -p results/tables results/figures
    Rscript ${projectDir}/scripts/07_external_validation_gse71729.R
    """
}

workflow {
    FETCH_AND_CURATE()
    SUBTYPING_AND_DESEQ2(FETCH_AND_CURATE.out.counts, FETCH_AND_CURATE.out.clinical)
    PATHWAY_AND_SURVIVAL(FETCH_AND_CURATE.out.counts, SUBTYPING_AND_DESEQ2.out.clinical_subtyped, SUBTYPING_AND_DESEQ2.out.deg_table)
    TME_DECONVOLUTION(FETCH_AND_CURATE.out.counts, SUBTYPING_AND_DESEQ2.out.clinical_subtyped)
    DRUG_REPURPOSING(SUBTYPING_AND_DESEQ2.out.deg_table)
    STROMAL_ADJUSTED_DESEQ2_AND_DRUGS(FETCH_AND_CURATE.out.counts, SUBTYPING_AND_DESEQ2.out.clinical_subtyped, SUBTYPING_AND_DESEQ2.out.deg_table, TME_DECONVOLUTION.out.est_scores)
    EXTERNAL_VALIDATION(STROMAL_ADJUSTED_DESEQ2_AND_DRUGS.out.deg_adj)
}