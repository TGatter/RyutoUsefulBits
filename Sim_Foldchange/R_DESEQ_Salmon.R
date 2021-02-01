#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

folder = args[1]
count = args[2]

library(DESeq2)
library(tximportData)
library(tximport)
library(readr)


setwd(paste(folder,"/quants", sep=""))
files_list = list.files()
files <- file.path(paste(folder,"/quants", sep=""), "", files_list, "quant.sf")

setwd(folder)

names(files) <- c(paste("CTR",seq(1,count,1), sep=""), paste("FC",seq(1,count,1), sep=""))

txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
condition = factor(c(rep("CTR", count), rep("FC", count)))
ExpDesign <- data.frame(row.names=colnames(txi.tx$counts), condition = condition)

dds <- DESeqDataSetFromTximport(txi.tx, ExpDesign, ~condition)
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","CTR","FC"))
resOrdered <- res[order(res$pvalue),]

resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig),  file="CTR_vs_FC_Salmon_0.1.csv")

resSig <- subset(resOrdered, padj < 0.05)
write.csv(as.data.frame(resSig),  file="CTR_vs_FC_Salmon_0.05.csv")

resSig <- subset(resOrdered, padj < 0.01)
write.csv(as.data.frame(resSig),  file="CTR_vs_FC_Salmon_0.01.csv")


