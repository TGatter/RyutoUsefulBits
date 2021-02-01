#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

folder =  args[1]
count =  args[2]

library(DESeq2)
library(tximportData)
library(tximport)
library(readr)


setwd(folder)

cts <- as.matrix(read.csv( "transcripts.count" ,sep="\t",row.names="Transcript"))
rn <- rownames(cts)
cts <- cts[,-1]
cts <- cts[,-1]
cts <- apply(cts, 2, as.integer)

cts[is.na(cts)] <- 0

rownames(cts) <- rn

coldata <- data.frame(condition=c(rep("CTR", count), rep("FC", count)))

rownames(coldata) <- c(paste("CTR",seq(1,count,1), sep=""), paste("FC",seq(1,count,1), sep=""))

condition = factor(c(rep("CTR", count), rep("FC", count)))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","CTR","FC"))
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig),  file="CTR_vs_FC_Ryuto_0.1.csv")

resSig <- subset(resOrdered, padj < 0.05)
write.csv(as.data.frame(resSig),  file="CTR_vs_FC_Ryuto_0.05.csv")

resSig <- subset(resOrdered, padj < 0.01)
write.csv(as.data.frame(resSig),  file="CTR_vs_FC_Ryuto_0.01.csv")
