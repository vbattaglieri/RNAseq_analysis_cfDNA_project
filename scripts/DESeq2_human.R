
##Analisi di HCT116 WT and DKO per cfDNA (03.05.24)


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


# carico sample_data con il design sperimentale per deseq2

sample_data <- read.table("samples_info.txt", header=T, sep="\t")
rownames(sample_data) <- sample_data$pop
sample_data$pop <- as.factor(sample_data$pop)
sample_data$run <- as.factor(sample_data$run)
sample_data$condition <- as.factor(sample_data$condition)

# carico i risultati di RSEM con le conte

files_DKO <- list.files("rsem", "*genes.results$", full.names = T)

# carico i nomi dei file di RSEM 

files_name_DKO <- list.files("rsem", "*genes.results$")

# rimuovo parte finale dei nomi

files_name_DKO <- str_remove(files_name_DKO, ".genes.results")

# aggiungo rownames

names(files_DDT) <- files_name_DDT

# importo le conte con txi (sia DKO che WT)

txi_DKO <- tximport(files_DKO, type = "rsem", txIn = TRUE, txOut = TRUE)

# modifico length dei numeri

txi_DKO$length[txi_DKO$length == 0] = 0.01

# Creo oggetto DEseq con sample_data infos e design con condition

dds_DKO <- DESeqDataSetFromTximport(txi_DKO, colData = sample_data, design = ~ condition)

# Relevel

dds_DKO$condition <- relevel(dds_DKO$condition, ref = "WT")

# Performo Deseq2

dds_DKO <- DESeq(dds_DKO)

# Performo LFCshrinkage

dds_DKO_lfcShrink <- lfcShrink(dds_DKO, coef="condition_DKO_vs_WT", type="apeglm")

# Stampo come tabella prima risultato di geni differenzialmente espressi

write.table(dds_DKO_lfcShrink, "dds_DKO_lfcShrink.txt", sep="\t", quote=FALSE, row.names=TRUE)

# Tolgo il . dopo ENSG 

rownames(dds_DKO_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_DKO_lfcShrink))

# Si crea la lista di annotazione nome del gene (symbol) e ensembl id 

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=rownames(dds_DKO_lfcShrink),
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# Aggiungo una colonna chiamata ENSEMBL 

ds_DKO_lfcShrink$ENSEMBL <- rownames(dds_DKO_lfcShrink)

#


