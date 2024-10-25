#!/usr/bin/env Rscript

suppressMessages(library("optparse"))

### Script to extract expression values of specific control genes 

option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix was originally produced using featureCounts or not. If true, rpkm values will be considered instead of cpm. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] genes_list metadata
                     'genes_list' is the path to a tab-delimited file containing the list of control genes to evaluate. Each gene must be reported in a row, with no header.
                     'metadata' is the path to a tab delimited file containing sample metadata. Sample names must be reported in the first column.",
                     option_list = option_list, prog = "control_genes",
                     description = "Extract expression values (rpkm or cpm) of specific control genes."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=2) {
  stop("Two arguments must be supplied (genes_list and metadata)", call.=FALSE)
}

genes_list = arguments$args[1]
metadata = arguments$args[2]
feature_counts = opt$featureCounts


# Importing expression values
if(feature_counts){
  data = read.delim("log.norm.rpkm.txt", h=T, row.names=1)
} else {
  data = read.delim("cpm.txt", h=T, row.names=1)
}

# Importing metadata and matching names
meta = read.delim(metadata, h=T, row.names=1, check.names=F)
meta = meta[rownames(meta)%in%colnames(data),]
meta = meta[match(colnames(data), rownames(meta)),]

# Importing control genes
control_genes = as.character(read.delim(genes_list, h=F)$V1)
control_genes = control_genes[control_genes%in%rownames(data)]
if(length(control_genes)==0){
  stop("No control genes were retrieved in the expression matrix. Check consistency between the expression data and the supplied list of genes.")
}

# Filtering data table
data_ctrl = data[control_genes,]
data_ctrl = t(data_ctrl)
data_ctrl = cbind(data_ctrl, meta)
data_ctrl = cbind(rownames(data_ctrl), data_ctrl)
colnames(data_ctrl)[1] = "sample"

write.table(data_ctrl, "control_genes_exprs.txt", quote=F, row.names=F, col.names=T, sep="\t")




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()

