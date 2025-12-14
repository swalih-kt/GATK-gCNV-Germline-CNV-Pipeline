#!/bin/bash
# ============================================================
# GATK gCNV Germline CNV Pipeline
# Author: Mohammed Swalih KT
# Reference: hg38 (GATK Resource Bundle)
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# 0. INPUT FILES & DIRECTORIES
# ------------------------------------------------------------

# Directory containing BAM files (sorted + indexed)
BAM_DIR="/path/to/bam_files"

# BED file with CNV regions of interest
BED_FILE="interval_chr_list_with_chr_updated.bed"

# Reference genome
REF="/home/mercy/Desktop/swalih_work/treesa/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

# Output directory
OUTDIR="gcnv_results"
mkdir -p ${OUTDIR}

# ------------------------------------------------------------
# 1. SORT AND INDEX BAM FILES (IF NOT DONE)
# ------------------------------------------------------------

# Sort BAM files
for bam in ${BAM_DIR}/*.bam; do
    samtools sort -o ${bam%.bam}.sorted.bam ${bam}
    samtools index ${bam%.bam}.sorted.bam
done

# ------------------------------------------------------------
# 2. PREPROCESS INTERVALS
# ------------------------------------------------------------

# 2.1 Preprocess intervals using BED file
gatk --java-options "-Xmx8G" PreprocessIntervals \
    -R ${REF} \
    -L ${BED_FILE} \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTDIR}/targets.preprocessed.interval_list

# 2.2 Genome-wide preprocessing (used for GC annotation)
gatk --java-options "-Xmx8G" PreprocessIntervals \
    -R ${REF} \
    --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTDIR}/grch38.preprocessed.interval_list

# ------------------------------------------------------------
# 3. COLLECT READ COUNTS (RUN FOR EACH BAM)
# ------------------------------------------------------------

mkdir -p ${OUTDIR}/counts

for bam in ${BAM_DIR}/*.sorted.bam; do
    sample=$(basename ${bam%.sorted.bam})

    gatk --java-options "-Xmx8G" CollectReadCounts \
        -L ${OUTDIR}/targets.preprocessed.interval_list \
        -R ${REF} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -I ${bam} \
        --format TSV \
        -O ${OUTDIR}/counts/${sample}.tsv
done

# ------------------------------------------------------------
# 4. ANNOTATE INTERVALS (GC CONTENT)
# ------------------------------------------------------------

gatk --java-options "-Xmx8G" AnnotateIntervals \
    -L ${OUTDIR}/grch38.preprocessed.interval_list \
    -R ${REF} \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTDIR}/annotated_intervals.tsv

# ------------------------------------------------------------
# 5. FILTER INTERVALS (GC BIAS & EXTREME COUNTS)
# ------------------------------------------------------------

gatk --java-options "-Xmx8G" FilterIntervals \
    -L ${OUTDIR}/grch38.preprocessed.interval_list \
    --annotated-intervals ${OUTDIR}/annotated_intervals.tsv \
    -I ${OUTDIR}/counts/*.tsv \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTDIR}/cohort.gc.filtered.interval_list

# ------------------------------------------------------------
# 6. DETERMINE CONTIG PLOIDY (COHORT MODE)
# ------------------------------------------------------------

gatk --java-options "-Xmx8G" DetermineGermlineContigPloidy \
    -L ${OUTDIR}/cohort.gc.filtered.interval_list \
    -I ${OUTDIR}/counts/*.tsv \
    --contig-ploidy-priors ploidy_priors.tsv \
    --output ${OUTDIR}/ploidy \
    --output-prefix cohort_ploidy \
    --verbosity DEBUG

# ------------------------------------------------------------
# 7. SCATTER INTERVALS FOR CNV CALLING
# ------------------------------------------------------------

mkdir -p ${OUTDIR}/scatter

gatk IntervalListTools \
    --INPUT ${OUTDIR}/cohort.gc.filtered.interval_list \
    --SUBDIVISION_MODE INTERVAL_COUNT \
    --SCATTER_CONTENT 15000 \
    --OUTPUT ${OUTDIR}/scatter/scatter_intervals

# ------------------------------------------------------------
# 8. GERMLINE CNV CALLING (COHORT MODE)
# ------------------------------------------------------------

gatk --java-options "-Xmx8G" GermlineCNVCaller \
    --run-mode COHORT \
    -L ${OUTDIR}/scatter/scatter_intervals.interval_list \
    -I ${OUTDIR}/counts/*.tsv \
    --contig-ploidy-calls ${OUTDIR}/ploidy \
    --annotated-intervals ${OUTDIR}/annotated_intervals.tsv \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output ${OUTDIR}/cohort_cnv \
    --output-prefix cohort_cnv \
    --verbosity DEBUG

# ------------------------------------------------------------
# 9. POSTPROCESS CNV CALLS (PER SAMPLE)
# ------------------------------------------------------------

# NOTE:
# sample-index = order of BAM files (0-based)
# Example: sample-index 1 → second BAM file

gatk --java-options "-Xmx8G" PostprocessGermlineCNVCalls \
    --model-shard-path ${OUTDIR}/cohort_cnv-model \
    --calls-shard-path ${OUTDIR}/cohort_cnv-calls \
    --allosomal-contig chrX \
    --allosomal-contig chrY \
    --contig-ploidy-calls ${OUTDIR}/ploidy \
    --sample-index 1 \
    --output-denoised-copy-ratios ${OUTDIR}/denoised_copy_ratios_sample1.vcf.gz \
    --output-genotyped-intervals ${OUTDIR}/genotyped_intervals_sample1.vcf.gz \
    --output-genotyped-segments ${OUTDIR}/genotyped_segments_sample1.vcf.gz \
    --sequence-dictionary ${REF%.fasta}.dict
