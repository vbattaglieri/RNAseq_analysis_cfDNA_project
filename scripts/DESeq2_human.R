
##Analysis of HCT116 WT and DKO for cfDNA release


#!/usr/bin/env Rscript

#load required library 

library(data.table)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(stringr)
library(pheatmap)
library(tximport)
library(tidyverse)
library(org.Hs.eg.db)


# load sample data containing the experimental design for DESeq2

sample_data <- read.table("samples_info.txt", header=T, sep="\t")
rownames(sample_data) <- sample_data$pop
sample_data$pop <- as.factor(sample_data$pop)
sample_data$run <- as.factor(sample_data$run)
sample_data$condition <- as.factor(sample_data$condition)

# Load RSEM results and counts

files_DKO <- list.files("rsem", "*genes.results$", full.names = T)

# Load RSEM filenames

files_name_DKO <- list.files("rsem", "*genes.results$")

# Format samplename

files_name_DKO <- str_remove(files_name_DKO, ".genes.results")

# add rownames

names(files_DDT) <- files_name_DDT

# Import counts with txi (both WT and DKO)

txi_DKO <- tximport(files_DKO, type = "rsem", txIn = TRUE, txOut = TRUE)

# modify number length

txi_DKO$length[txi_DKO$length == 0] = 0.01

# Create DEseq object with sample_data info and design with condition

dds_DKO <- DESeqDataSetFromTximport(txi_DKO, colData = sample_data, design = ~ condition)

# Relevel

dds_DKO$condition <- relevel(dds_DKO$condition, ref = "WT")

# Perform Deseq2

dds_DKO <- DESeq(dds_DKO)

# Perform LFCshrinkage

dds_DKO_lfcShrink <- lfcShrink(dds_DKO, coef="condition_DKO_vs_WT", type="apeglm")

# Print table

write.table(dds_DKO_lfcShrink, "dds_DKO_lfcShrink.txt", sep="\t", quote=FALSE, row.names=TRUE)

# Format ENSG names

rownames(dds_DKO_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_DKO_lfcShrink))

# An annotation list linking gene symbols to Ensembl IDs is created

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=rownames(dds_DKO_lfcShrink),
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# Add column to ENSEMBL 

ds_DKO_lfcShrink$ENSEMBL <- rownames(dds_DKO_lfcShrink)

#


