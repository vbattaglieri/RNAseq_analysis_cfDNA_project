#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=cfDNA_RNA_HCT116
#SBATCH --error=cfDNA_RNA_HCT116_error.log
#SBATCH --output=cfDNA_RNA_HCT116_output.log
#SBATCH --cpus-per-task=40
#SBATCH --mem=150G
#SBATCH --time=96:00:00


# Usage: sbatch launcher_RNA.sh
# Arg 1 to the pipeline = directory containing the paired-end *_1.fastq.gz / *_2.fastq.gz files (no trailing slash)
bash /path/to/scripts/rna_seq_pipeline_QC_rsem_GRCh38_gencode.v44.sh /path/to/rawdata/RNA/HCT116




