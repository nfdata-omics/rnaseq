#!/usr/bin/env Rscript

library(edgeR)
library(limma)

# Importing raw counts - from Salmon
counts = read.delim("gene_counts.tsv", h=T, row.names=1) # importing data using ENSG
counts$gene_name <- make.unique(counts$gene_name) # make unique gene_name
rownames(counts) = counts$gene_name # switch to gene_name
counts = counts[,-1]

y_all = DGEList(counts=counts)
cpm = cpm(y_all, log=FALSE)
write.table(cpm, "cpm.txt", quote=F, row.names=T, col.names=T, sep="\t")

size_factors = data.frame(samples=rownames(y_all$samples), lib.size=y_all$samples$lib.size)
write.table(size_factors, "lib_size_factors.cpm.txt", quote=F, row.names=T, col.names=T, sep="\t")
