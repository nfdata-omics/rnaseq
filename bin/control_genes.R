#!/usr/bin/env Rscript

suppressMessages(library("optparse"))

### Script to extract expression values of specific control genes 

option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix was originally produced using featureCounts or not. If true, rpkm values will be considered instead of cpm. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] genes_list
                     genes_list is the path to a tab-delimited file containing the list of control genes to evaluate. Each gene must be reported in a row, with no header.",
                     option_list = option_list, prog = "control_genes",
                     description = "Extract expression values (rpkm or cpm) of specific control genes."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=1) {
  stop("One argument must be supplied (control_genes)", call.=FALSE)
}

genes_list = arguments$args[1]
feature_counts = opt$featureCounts


# Importing expression values
if(feature_counts){
  data = read.delim("log.norm.rpkm.txt", h=T, row.names=1)
} else {
  data = read.delim("cpm.txt", h=T, row.names=1)
}

# Importing control genes
control_genes = as.character(read.delim(genes_list, h=F)$V1)
control_genes = control_genes[control_genes%in%rownames(data)]
if(length(control_genes)==0){
  stop("No control genes were retrieved in the expression matrix. Check consistency between the expression data and the supplied list of genes.")
}

# Filtering data table
data_ctrl = data[control_genes,]
data_ctrl = cbind(rownames(data_ctrl), data_ctrl)
colnames(data_ctrl)[1] = "gene_name"
write.table(data_ctrl, "control_genes_exprs.txt", quote=F, row.names=F, col.names=T, sep="\t")




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()

