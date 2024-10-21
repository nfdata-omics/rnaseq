#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("edgeR"))
suppressMessages(library("limma"))

### Script to perform raw counts normalization into cpm and/or rpkm

option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix has been produced using featureCounts or not. If true, gene names are supposed to be reported in the first column, and gene informations (Chr/Start/End/Strand/Length) are supposed to be reported in the 2-6 columns. If false, only the first column containing gene names is expected. If true, both cpm and rpkm will be calculated; if false, only cpm will be calculated. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] input_data",
                     option_list = option_list, prog = "count_norm",
                     description = "Perform raw counts normalization into cpm and/or rpkm. 
                     'input_data' is the path of a tab delimited file containing the raw counts matrix."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=1) {
  stop("One argument must be supplied (input_data)", call.=FALSE)
}

input_data = arguments$args[1]
feature_counts = opt$featureCounts


# Importing raw counts - TODO: decomment this line after Carninci
counts = read.delim(input_data, h=T, row.names=1, check.names=F)

if(feature_counts){

	y_all = DGEList(counts=counts[,6:ncol(counts)], genes=counts[,1:5])
	# Cpm
	cpm = cpm(y_all, log=FALSE)
	cpm = cbind(rownames(cpm), cpm)
	colnames(cpm)[1] = "gene_name"
	write.table(cpm, "cpm.txt", quote=F, row.names=F, col.names=T, sep="\t")
	# Log-norm-rpkm (TMM normalization)
	y_all = calcNormFactors(y_all, method="TMM")
	log_norm_rpkm = rpkm(y_all, log=T, gene.length=y_all$genes$Length)
	log_norm_rpkm = cbind(rownames(log_norm_rpkm), log_norm_rpkm)
        colnames(log_norm_rpkm)[1] = "gene_name"
	write.table(log_norm_rpkm, "log.norm.rpkm.txt", quote=F, row.names=F, col.names=T, sep="\t")
	# Lib size and normalization factor
	size_factors = data.frame(samples=rownames(y_all$samples), lib.size=y_all$samples$lib.size, norm.factors=y_all$samples$norm.factors)
	write.table(size_factors, "lib_size_factors.txt", quote=F, row.names=F, col.names=T, sep="\t")

} else {

	# TODO: This is a temporary custom script for Carninci project fantom6, then delete this part
	# Importing raw counts - from Salmon
	#counts$gene_name <- make.unique(counts$gene_name) # make unique gene_name
	#rownames(counts) = counts$gene_name # switch to gene_name
	#counts = counts[,-1]
	### \end of Carninci custom script

	y_all = DGEList(counts=counts)
	# Cpm
	cpm = cpm(y_all, log=FALSE)
	cpm = cbind(rownames(cpm), cpm)
        colnames(cpm)[1] = "gene_name"
	write.table(cpm, "cpm.txt", quote=F, row.names=F, col.names=T, sep="\t")
	# Lib size
	size_factors = data.frame(samples=rownames(y_all$samples), lib.size=y_all$samples$lib.size)
	write.table(size_factors, "lib_size_factors.txt", quote=F, row.names=F, col.names=T, sep="\t")

}





########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()
