# RNA-seq: HCT116 WT vs DKO

Analysis code for a bulk RNA-seq comparison of the colorectal cancer cell line
**HCT116 wild-type (WT)** versus its **DNMT1/DNMT3B double-knockout (DKO)**
derivative, 3 biological replicates per condition.

This repository showcases the **scripts** used in the study — from raw-read QC
and quantification through differential expression and gene-set enrichment.
It is **code-only**: raw sequencing data and large intermediate files are not
tracked here (see [Data availability](#data-availability)).

> Part of a published project. If you use or refer to this code, please cite the
> associated publication *(add citation / DOI here)*.

## Experimental design

| Sample            | Condition |
|-------------------|-----------|
| HCT116_WT_Rep1-3  | WT        |
| HCT116_DKO_Rep1-3 | DKO       |

Full design table: [`metadata/samples_info.txt`](metadata/samples_info.txt).
Comparison: **DKO vs WT** (WT as reference level).

## Pipeline overview

```
FASTQ ──▶ FastQC ──▶ RSEM (STAR) quantification ──▶ RSeQC ──▶ DESeq2 ──▶ fgsea
         (QC)        GRCh38 / GENCODE v44          (QC)      (DE)      (GSEA)
```

| Stage | Tool | Script |
|-------|------|--------|
| QC | FastQC / MultiQC | `scripts/rna_seq_pipeline_QC_rsem_GRCh38_gencode.v44.sh` |
| Quantification | RSEM + STAR (paired-end, reverse-stranded) | same as above |
| Alignment QC | RSeQC (`bam_stat`, `infer_experiment`, `geneBody_coverage`) | same as above |
| Job submission | SLURM | `scripts/launcher_RNA.sh` |
| Differential expression | DESeq2 (apeglm LFC shrinkage) | `scripts/DESeq2_human.R`, `scripts/DESeq2.R` |
| Gene-set enrichment | fgsea (MSigDB) | `scripts/DESeq2.R` |

**Reference:** human genome **GRCh38**, annotation **GENCODE v44**.

## Repository layout

```
.
├── scripts/       # QC + quantification pipeline (bash) and DE/GSEA analysis (R)
├── metadata/      # experimental design (samples_info.txt)
└── results/       # selected final figures (volcano, DE heatmap, GSEA)
```

## Requirements

- **Shell pipeline:** FastQC, MultiQC, RSEM, STAR, samtools, RSeQC
  (a conda env named `rnaseq` is assumed by the scripts).
- **R (≥ 4.3):** DESeq2, apeglm, tximport, fgsea, EnhancedVolcano, pheatmap,
  ggplot2, ggrepel, org.Hs.eg.db, tidyverse, data.table, stringr.

## Usage

1. **Configure paths.** Edit the `BASE_DIR` / `CONDA_ENV` block at the top of
   `scripts/rna_seq_pipeline_QC_rsem_GRCh38_gencode.v44.sh` and the FASTQ path in
   `scripts/launcher_RNA.sh` for your environment. Build the RSEM/STAR index for
   GRCh38 + GENCODE v44 beforehand.

2. **Run QC + quantification** (SLURM):
   ```bash
   sbatch scripts/launcher_RNA.sh
   ```
   This produces per-sample `rsem/*.genes.results`.

3. **Differential expression + GSEA** (from the analysis directory containing
   `rsem/` and `samples_info.txt`):
   ```bash
   Rscript scripts/DESeq2_human.R
   ```

## Data availability

Raw and processed sequencing data are deposited in a public repository
*(add GEO/SRA/ArrayExpress accession here)*. MSigDB gene-set (`.gmt`) files are
**not redistributed** — download them from
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) under their license.

## License

Code released under the [MIT License](LICENSE).
