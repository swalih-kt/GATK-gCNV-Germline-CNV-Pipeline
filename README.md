# GATK gCNV Germline CNV Analysis Pipeline (WES/WGS)

This repository documents a **GATK gCNV-based workflow for germline Copy Number Variant (CNV) detection** using BAM files.  
The pipeline follows **GATK Best Practices** and includes interval preprocessing, read count collection, GC filtering, contig ploidy estimation, CNV calling, and post-processing.

This workflow is suitable for **WES/WGS cohort-based CNV analysis** in human samples (hg38).

---

## 📌 Overview of the Pipeline

The pipeline performs the following major steps:

1. BAM file preprocessing (sorting & indexing)
2. Interval preparation (BED → interval lists)
3. Read count collection per interval
4. Interval annotation and GC filtering
5. Contig ploidy determination (autosomal + sex chromosomes)
6. Germline CNV calling (cohort mode)
7. Postprocessing and sample-level CNV extraction

---

## 🧬 Input Requirements

### 1. BAM Files
- All BAM files must be:
  - Coordinate-sorted
  - Indexed (`.bam.bai`)
  - Located in the **same directory**

Example:
```bash
sample1.bam
sample1.bam.bai
sample2.bam
sample2.bam.bai
