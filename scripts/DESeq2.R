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
library(org.Mm.eg.db)


#dataset info loading
sample_data <- read.table("results/samples_info.txt", header=T, sep="\t")
rownames(sample_data) <- sample_data$pop
sample_data$pop <- as.factor(sample_data$pop)
sample_data$run <- as.factor(sample_data$run)
sample_data$condition <- as.factor(sample_data$condition)

#sample_data

sample_data_DD <- read.table("results/samples_info_DD.txt", header=T, sep="\t")
rownames(sample_data_DD) <- sample_data_DD$pop
sample_data_DD$pop <- as.factor(sample_data_DD$pop)
sample_data_DD$condition <- as.factor(sample_data_DD$condition)

sample_data_DDT <- read.table("results/samples_info_DDT.txt", header=T, sep="\t")
rownames(sample_data_DDT) <- sample_data_DDT$pop
sample_data_DDT$pop <- as.factor(sample_data_DDT$pop)
sample_data_DDT$condition <- as.factor(sample_data_DDT$condition)

sample_data_TMZ <- read.table("results/samples_info_TMZ.txt", header=T, sep="\t")
rownames(sample_data_TMZ) <- sample_data_TMZ$pop
sample_data_TMZ$pop <- as.factor(sample_data_TMZ$pop)
sample_data_TMZ$condition <- as.factor(sample_data_TMZ$condition)

sample_data_FOXIR <- read.table("results/samples_info_FOXIR.txt", header=T, sep="\t")
rownames(sample_data_FOXIR) <- sample_data_FOXIR$pop
sample_data_FOXIR$pop <- as.factor(sample_data_FOXIR$pop)
sample_data_FOXIR$condition <- as.factor(sample_data_FOXIR$condition)

sample_data_TEMIRI <- read.table("results/samples_info_TEMIRI.txt", header=T, sep="\t")
rownames(sample_data_TEMIRI) <- sample_data_TEMIRI$pop
sample_data_TEMIRI$pop <- as.factor(sample_data_TEMIRI$pop)
sample_data_TEMIRI$condition <- as.factor(sample_data_TEMIRI$condition)


#genes_results

files_DDT <- list.files("results_DDT__NT", "*genes.results$", full.names = T)
files_DD <- list.files("results_DD__NT", "*genes.results$", full.names = T)
files_TMZ <- list.files("results_TMZ__NT", "*genes.results$", full.names = T)
files_FOXIR <- list.files("results_FOXIR__NT", "*genes.results$", full.names = T)
files_TEMIRI <- list.files("results_TEMIRI__NT", "*genes.results$", full.names = T)

#format

files_name_DDT <- list.files("results_DDT__NT", "*genes.results$")
files_name_DD <- list.files("results_DD__NT", "*genes.results$")
files_name_TMZ <- list.files("results_TMZ__NT", "*genes.results$")
files_name_FOXIR <- list.files("results_FOXIR__NT", "*genes.results$")
files_name_TEMIRI <- list.files("results_TEMIRI__NT", "*genes.results$")


#files_name <- list.files("results", "*genes.results$")
files_name_DDT <- str_remove(files_name_DDT, ".genes.results")
files_name_DD <- str_remove(files_name_DD, ".genes.results")
files_name_TMZ <- str_remove(files_name_TMZ, ".genes.results")
files_name_FOXIR <- str_remove(files_name_FOXIR, ".genes.results")
files_name_TEMIRI <- str_remove(files_name_TEMIRI, ".genes.results")

#import
names(files_DDT) <- files_name_DDT
names(files_DD) <- files_name_DD
names(files_TMZ) <- files_name_TMZ
names(files_FOXIR) <- files_name_FOXIR
names(files_TEMIRI) <- files_name_TEMIRI

txi_DDT <- tximport(files_DDT, type = "rsem", txIn = TRUE, txOut = TRUE)
txi_DD <- tximport(files_DD, type = "rsem", txIn = TRUE, txOut = TRUE)
txi_TMZ <- tximport(files_TMZ, type = "rsem", txIn = TRUE, txOut = TRUE)
txi_FOXIR <- tximport(files_FOXIR, type = "rsem", txIn = TRUE, txOut = TRUE)
txi_TEMIRI <- tximport(files_TEMIRI, type = "rsem", txIn = TRUE, txOut = TRUE)

txi_DDT$length[txi_DDT$length == 0] = 0.01
txi_DD$length[txi_DD$length == 0] = 0.01
txi_TMZ$length[txi_TMZ$length == 0] = 0.01
txi_FOXIR$length[txi_FOXIR$length == 0] = 0.01
txi_TEMIRI$length[txi_TEMIRI$length == 0] = 0.01


#countData_DD <- txi_rsem$counts[,c(1,2,3,10,11,12)]
#countData_DDT <- txi_rsem$counts[,c(4,5,6,10,11,12)]
#countData_TMZ <- txi_rsem$counts[,c(16,17,18,10,11,12)]
#countData_FOXIR <- txi_rsem$counts[,c(7,8,9,10,11,12)]
#countData_TEMIRI <- txi_rsem$counts[,c(13,14,15,10,11,12)]

#Independent DESeq

dds_DD <- DESeqDataSetFromTximport(txi_DD, colData = sample_data_DD, design = ~ condition)
dds_DDT <- DESeqDataSetFromTximport(txi_DDT, colData = sample_data_DDT, design = ~ condition)
dds_TMZ <- DESeqDataSetFromTximport(txi_TMZ, colData = sample_data_TMZ, design = ~ condition)
dds_FOXIR <- DESeqDataSetFromTximport(txi_FOXIR, colData = sample_data_FOXIR, design = ~ condition)
dds_TEMIRI <- DESeqDataSetFromTximport(txi_TEMIRI, colData = sample_data_TEMIRI, design = ~ condition)

dds_DD$condition <- relevel(dds_DD$condition, ref = "NT")
dds_DDT$condition <- relevel(dds_DDT$condition, ref = "NT")
dds_TMZ$condition <- relevel(dds_TMZ$condition, ref = "NT")
dds_FOXIR$condition <- relevel(dds_FOXIR$condition, ref = "NT")
dds_TEMIRI$condition <- relevel(dds_TEMIRI$condition, ref = "NT")

dds_DD <- DESeq(dds_DD)
dds_DDT <- DESeq(dds_DDT)
dds_TMZ <- DESeq(dds_TMZ)
dds_FOXIR <- DESeq(dds_FOXIR)
dds_TEMIRI <- DESeq(dds_TEMIRI)

resultsNames(dds_DD)
resultsNames(dds_DDT)
resultsNames(dds_TMZ)
resultsNames(dds_FOXIR)
resultsNames(dds_TEMIRI)

dds_DD_lfcShrink <- lfcShrink(dds_DD, coef="condition_DD_vs_NT", type="apeglm")
dds_DDT_lfcShrink <- lfcShrink(dds_DDT, coef="condition_DDT_vs_NT", type="apeglm")
dds_TMZ_lfcShrink <- lfcShrink(dds_TMZ, coef="condition_TMZ_vs_NT", type="apeglm")
dds_FOXIR_lfcShrink <- lfcShrink(dds_FOXIR, coef="condition_FOXIR_vs_NT", type="apeglm")
dds_TEMIRI_lfcShrink <- lfcShrink(dds_TEMIRI, coef="condition_TEMIRI_vs_NT", type="apeglm")

#Save the differential expression results for each condition compared with T0
#results_DDT_vs_NT <- results(dds_DDT_lfcShrink, contrast = c("condition", "DDT", "NT"))
#results_FOXIR_vs_NT <- results(dds_FOXIR_lfcShrink, contrast = c("condition", "FOXIR", "NT"))
#results_DD_vs_NT <- results(dds_DD_lfcShrink, contrast = c("condition", "DD", "NT"))
#results_TMZ_vs_NT <- results(dds_TMZ_lfcShrink, contrast = c("condition", "TMZ", "NT"))
#results_TEMIRI_vs_NT <- results(dds_TEMIRI_lfcShrink, contrast = c("condition", "TEMIRI", "NT"))

#write.table(results_DDT_vs_NT, "res_dds_chemo_priming_DD_DESeq.txt", sep="\t", quote=FALSE, row.names=TRUE)
#write.table(results_FOXIR_vs_NT, "res_dds_chemo_priming_DDT_DESeq.txt", sep="\t", quote=FALSE, row.names=TRUE)
#write.table(results_DD_vs_NT, "res_dds_chemo_priming_FOXIR_DESeq.txt", sep="\t", quote=FALSE, row.names=TRUE)
#write.table(results_TMZ_vs_NT, "res_dds_chemo_priming_TEMIRI_DESeq.txt", sep="\t", quote=FALSE, row.names=TRUE)
#write.table(results_TEMIRI_vs_NT, "res_dds_chemo_priming_TMZ_DESeq.txt", sep="\t", quote=FALSE, row.names=TRUE)


rownames(dds_DD_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_DD_lfcShrink))
rownames(dds_DDT_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_DDT_lfcShrink))
rownames(dds_TMZ_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_TMZ_lfcShrink))
rownames(dds_FOXIR_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_FOXIR_lfcShrink))
rownames(dds_TEMIRI_lfcShrink) <-  gsub("\\.[0-9]*$", "", rownames(dds_TEMIRI_lfcShrink))

length(rownames(dds_DD_lfcShrink))
length(rownames(dds_DDT_lfcShrink))
length(rownames(dds_TMZ_lfcShrink))
length(rownames(dds_FOXIR_lfcShrink))
length(rownames(dds_TEMIRI_lfcShrink))

ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
                                    key=rownames(dds_DD_lfcShrink), 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
#ens2symbol <- as_tibble(ens2symbol)

#ens2symbol

dds_DD_lfcShrink$ENSEMBL <- rownames(dds_DD_lfcShrink)
dds_DDT_lfcShrink$ENSEMBL <- rownames(dds_DDT_lfcShrink)
dds_TMZ_lfcShrink$ENSEMBL <- rownames(dds_TMZ_lfcShrink)
dds_FOXIR_lfcShrink$ENSEMBL <- rownames(dds_FOXIR_lfcShrink)
dds_TEMIRI_lfcShrink$ENSEMBL <- rownames(dds_TEMIRI_lfcShrink)

dds_DD_lfcShrink <- as.data.frame(dds_DD_lfcShrink)
dds_DDT_lfcShrink <- as.data.frame(dds_DDT_lfcShrink)
dds_TMZ_lfcShrink <- as.data.frame(dds_TMZ_lfcShrink)
dds_FOXIR_lfcShrink <- as.data.frame(dds_FOXIR_lfcShrink)
dds_TEMIRI_lfcShrink <- as.data.frame(dds_TEMIRI_lfcShrink)


dds_DD_lfcShrink_merge <- merge(dds_DD_lfcShrink, ens2symbol, by="ENSEMBL")
dds_DDT_lfcShrink_merge <- merge(dds_DDT_lfcShrink, ens2symbol, by="ENSEMBL")
dds_TMZ_lfcShrink_merge <- merge(dds_TMZ_lfcShrink, ens2symbol, by="ENSEMBL")
dds_FOXIR_lfcShrink_merge <- merge(dds_FOXIR_lfcShrink, ens2symbol, by="ENSEMBL")
dds_TEMIRI_lfcShrink_merge <- merge(dds_TEMIRI_lfcShrink, ens2symbol, by="ENSEMBL")


dds_DD_lfcShrink_merge <- na.omit(dds_DD_lfcShrink_merge[, c("SYMBOL", "log2FoldChange")])
dds_DDT_lfcShrink_merge <- na.omit(dds_DDT_lfcShrink_merge[, c("SYMBOL", "log2FoldChange")])
dds_TMZ_lfcShrink_merge <- na.omit(dds_TMZ_lfcShrink_merge[, c("SYMBOL", "log2FoldChange")])
dds_FOXIR_lfcShrink_merge <- na.omit(dds_FOXIR_lfcShrink_merge[, c("SYMBOL", "log2FoldChange")])
dds_TEMIRI_lfcShrink_merge <- na.omit(dds_TEMIRI_lfcShrink_merge[, c("SYMBOL", "log2FoldChange")])

dds_DD_lfcShrink_merge[!duplicated(dds_DD_lfcShrink_merge),]
dds_DDT_lfcShrink_merge[!duplicated(dds_DDT_lfcShrink_merge),]
dds_TMZ_lfcShrink_merge[!duplicated(dds_TMZ_lfcShrink_merge),]
dds_FOXIR_lfcShrink_merge[!duplicated(dds_FOXIR_lfcShrink_merge),]
dds_TEMIRI_lfcShrink_merge[!duplicated(dds_TEMIRI_lfcShrink_merge),]

dds_DD_lfcShrink_array <- dds_DD_lfcShrink_merge$log2FoldChange
dds_DDT_lfcShrink_array <- dds_DDT_lfcShrink_merge$log2FoldChange
dds_TMZ_lfcShrink_array <- dds_TMZ_lfcShrink_merge$log2FoldChange
dds_FOXIR_lfcShrink_array <- dds_FOXIR_lfcShrink_merge$log2FoldChange
dds_TEMIRI_lfcShrink_array <- dds_TEMIRI_lfcShrink_merge$log2FoldChange

names(dds_DD_lfcShrink_array) <- dds_DD_lfcShrink_merge$SYMBOL
names(dds_DDT_lfcShrink_array) <- dds_DDT_lfcShrink_merge$SYMBOL
names(dds_TMZ_lfcShrink_array) <- dds_TMZ_lfcShrink_merge$SYMBOL
names(dds_FOXIR_lfcShrink_array) <- dds_FOXIR_lfcShrink_merge$SYMBOL
names(dds_TEMIRI_lfcShrink_array) <- dds_TEMIRI_lfcShrink_merge$SYMBOL


#res_dds_DD_lfcShrink <- inner_join(dds_DD_lfcShrink, ens2symbol, by=c("row"="ENSEMBL"))
#res_dds_DDT_lfcShrink <- inner_join(dds_DDT_lfcShrink, ens2symbol, by=c("row"="ENSEMBL"))
#res_dds_TMZ_lfcShrink <- inner_join(dds_TMZ_lfcShrink, ens2symbol, by=c("row"="ENSEMBL"))
#res_dds_FOXIR_lfcShrink <- inner_join(dds_FOXIR_lfcShrink, ens2symbol, by=c("row"="ENSEMBL"))
#res_dds_TEMIRI_lfcShrink <- inner_join(dds_TEMIRI_lfcShrink, ens2symbol, by=c("row"="ENSEMBL"))

library(fgsea)

pathways.hallmark <- gmtPathways("mh.all.v2023.2.Mm.symbols.gmt")

fgseaRes_DD <- fgsea(pathways=pathways.hallmark, stats=dds_DD_lfcShrink_array, nperm=1000)
fgseaRes_DDT <- fgsea(pathways=pathways.hallmark, stats=dds_DDT_lfcShrink_array, nperm=1000)
fgseaRes_TMZ <- fgsea(pathways=pathways.hallmark, stats=dds_TMZ_lfcShrink_array, nperm=1000)
fgseaRes_TEMIRI <- fgsea(pathways=pathways.hallmark, stats=dds_FOXIR_lfcShrink_array, nperm=1000)
fgseaRes_FOXIR <- fgsea(pathways=pathways.hallmark, stats=dds_TEMIRI_lfcShrink_array, nperm=1000)

fgseaResTidy_DD <- fgseaRes_DD %>%
	as_tibble() %>%                                     
	arrange(desc(NES))

fgseaResTidy_DDT <- fgseaRes_DDT %>%
        as_tibble() %>%
        arrange(desc(NES))

fgseaResTidy_TMZ <- fgseaRes_TMZ %>%
        as_tibble() %>%
        arrange(desc(NES))

fgseaResTidy_TEMIRI <- fgseaRes_TEMIRI %>%
        as_tibble() %>%
        arrange(desc(NES))

fgseaResTidy_FOXIR <- fgseaRes_FOXIR %>%
        as_tibble() %>%
        arrange(desc(NES))




pdf("fgseaResTidy_DD.pdf", width=10)
ggplot(fgseaResTidy_DD, aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill=padj<0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA DD") +
	theme_minimal()
dev.off()


pdf("fgseaResTidy_DDT.pdf", width=10)
ggplot(fgseaResTidy_DDT, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA DDT") +
        theme_minimal()
dev.off()

pdf("fgseaResTidy_TMZ.pdf", width=10)
ggplot(fgseaResTidy_TMZ, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA TMZ") +
        theme_minimal()
dev.off()

pdf("fgseaResTidy_TEMIRI.pdf", width=10)
ggplot(fgseaResTidy_TEMIRI, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA TEMIRI") +
        theme_minimal()
dev.off()

pdf("fgseaResTidy_FOXIR.pdf", width=10)
ggplot(fgseaResTidy_FOXIR, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA FOXIR") +
        theme_minimal()
dev.off()


#m2.cp.reactome.v2023.2.Mm.symbols.gmt


#pathways.reactome <- gmtPathways("m2.cp.reactome.v2023.2.Mm.symbols.gmt")

#fgseaRes_DD <- fgsea(pathways=pathways.reactome, stats=dds_DD_lfcShrink_array, nperm=1000)
#fgseaRes_DDT <- fgsea(pathways=pathways.reactome, stats=dds_DDT_lfcShrink_array, nperm=1000)
#fgseaRes_TMZ <- fgsea(pathways=pathways.reactome, stats=dds_TMZ_lfcShrink_array, nperm=1000)
#fgseaRes_TEMIRI <- fgsea(pathways=pathways.reactome, stats=dds_FOXIR_lfcShrink_array, nperm=1000)
#fgseaRes_FOXIR <- fgsea(pathways=pathways.reactome, stats=dds_TEMIRI_lfcShrink_array, nperm=1000)

#fgseaResTidy_DD <- fgseaRes_DD %>%
#        as_tibble() %>%
#        arrange(desc(NES))

#fgseaResTidy_DDT <- fgseaRes_DDT %>%
#        as_tibble() %>%
#        arrange(desc(NES))

#fgseaResTidy_TMZ <- fgseaRes_TMZ %>%
#        as_tibble() %>%
#        arrange(desc(NES))

#fgseaResTidy_TEMIRI <- fgseaRes_TEMIRI %>%
#        as_tibble() %>%
#        arrange(desc(NES))

#fgseaResTidy_FOXIR <- fgseaRes_FOXIR %>%
#        as_tibble() %>%
#        arrange(desc(NES))


#pdf("fgseaResTidy_DD_Reactome.pdf", width=20, height=80)
#ggplot(fgseaResTidy_DD, aes(reorder(pathway, NES), NES)) +
#        geom_col(aes(fill=padj<0.05)) +
#        coord_flip() +
#        labs(x="Pathway", y="Normalized Enrichment Score", title="Reactome NES from GSEA DD") +
#        theme_minimal()
#dev.off()


#pdf("fgseaResTidy_DDT_Reactome.pdf", width=20, height=80)
#ggplot(fgseaResTidy_DDT, aes(reorder(pathway, NES), NES)) +
#        geom_col(aes(fill=padj<0.05)) +
#        coord_flip() +
#        labs(x="Pathway", y="Normalized Enrichment Score", title="Reactome NES from GSEA DDT") +
#        theme_minimal()
#dev.off()

#pdf("fgseaResTidy_TMZ_Reactome.pdf", width=20, height=80)
#ggplot(fgseaResTidy_TMZ, aes(reorder(pathway, NES), NES)) +
#        geom_col(aes(fill=padj<0.05)) +
#        coord_flip() +
#        labs(x="Pathway", y="Normalized Enrichment Score", title="Reactome NES from GSEA TMZ") +
#        theme_minimal()
#dev.off()

#pdf("fgseaResTidy_TEMIRI_Reactome.pdf", width=20, height=80)
#ggplot(fgseaResTidy_TEMIRI, aes(reorder(pathway, NES), NES)) +
#        geom_col(aes(fill=padj<0.05)) +
#        coord_flip() +
#        labs(x="Pathway", y="Normalized Enrichment Score", title="Reactome NES from GSEA TEMIRI") +
#        theme_minimal()
#dev.off()

#pdf("fgseaResTidy_FOXIR_Reactome.pdf", width=20, height=80)
#ggplot(fgseaResTidy_FOXIR, aes(reorder(pathway, NES), NES)) +
#        geom_col(aes(fill=padj<0.05)) +
#        coord_flip() +
#        labs(x="Pathway", y="Normalized Enrichment Score", title="Reactome NES from GSEA FOXIR") +
#        theme_minimal()
#dev.off()

print("End of analysis")

q()






