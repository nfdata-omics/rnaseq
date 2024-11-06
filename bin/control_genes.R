#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("pheatmap"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))

### Script to extract expression values of specific control genes 

option_list <- list(
  #make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix was originally produced using featureCounts or not. If true, rpkm values will be considered instead of cpm. [default \"%default\"]"),
  #I think it could be better to explicitly specify the name of the input table
  make_option(c("-a", "--annotations"), action="store", type="character", default="", help="The name(s) of categorical variables to highlight on the top of the heatmap (as annotation column bars). They must be the names of the corresponding metadata columns. If multiple names are provided, they must be comma-separated with no blank spaces (e.g. genotype,treatment). [default \"%default\"]"),
  make_option(c("-d", "--dendrogram"), action="store", type="integer", default=3, help="Whether to perform clustering and show the corresponding dendrogram on rows/genes ('1'), on columns/samples ('2'), both ('3') or none ('0'). [default \"%default\"]"),
  make_option(c("-s", "--scale_rows"), action="store_true", default=FALSE, help="Whether the values should be centered and scaled in the row direction. If true, genes with constant expression will be removed. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] input_data genes_list metadata
                     'input_data' is the path to a tab-delimited file containing the expression values (ideally cpm or rpkm if possible) from which the expression levels of selected control genes are extracted. Samples must be reported on different columns, genes on different rows.
                     'genes_list' is the path to a tab-delimited file containing the list of control genes to evaluate. Each gene must be reported in a different row, with no header.
                     'metadata' is the path to a tab delimited file containing sample metadata. Sample names must be reported in the first column.",
                     option_list = option_list, prog = "control_genes",
                     description = "Extract the expression values of specific control genes, and produce the corresponding heatmap."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=3) {
  stop("Three arguments must be supplied (input_data, genes_list and metadata)", call.=FALSE)
}

input_data = arguments$args[1]
genes_list = arguments$args[2]
metadata = arguments$args[3]
#feature_counts = opt$featureCounts
annot = str_split_1(opt$annotations, ",")
dend = opt$dendrogram
scale_rows = opt$scale_rows


# Importing expression values
# if(feature_counts){
#   data = read.delim("log.norm.rpkm.txt", h=T, row.names=1)
# } else {
#   data = read.delim("cpm.txt", h=T, row.names=1)
# }
data = read.delim(input_data, h=T, row.names=1)

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

# Excluding constant genes if scaling on row
if(scale_rows){
  var = apply(data_ctrl, 1, var)
  data_ctrl = data_ctrl[var!=0,]
  scale_value = "row"
} else{
  scale_value = "none"
}

# Annotations for heatmap
if(annot[1]!=""){
  if(sum(!annot%in%colnames(meta))>0){
    stop("Annotations for categorical variables must be names of metadata columns. If more than one variable are provided, they must be comma-separated with no blank spaces (e.g. treatment,group)")
  }
  col_annot = as.data.frame(meta[,annot])
  colnames(col_annot) = annot
  rownames(col_annot) = rownames(meta)
  for(i in 1:ncol(col_annot)){
    col_annot[,i] = as.factor(col_annot[,i])
  }
}


# Handling clustering options
if(dend==0) {clRows=F; clCols=F}
if(dend==1) {clRows=T; clCols=F}
if(dend==2) {clRows=F; clCols=T}
if(dend==3) {clRows=T; clCols=T}
if(dend>3) {stop("The dendrogram option must be set to an integer number ranging from 0 to 3.")}

# Plot
pdf("heatmap_control_genes.pdf", width=8, height=10)
if(annot[1]!=""){
  pheatmap(data_ctrl, main="Heatmap of selected genes", display_numbers=F, scale=scale_value, cluster_rows=clRows, cluster_cols=clCols, annotation_col=col_annot)
} else {
  pheatmap(data_ctrl, main="Heatmap of selected genes", display_numbers=F, scale=scale_value, cluster_rows=clRows, cluster_cols=clCols)
}
dev.off()


# Saving data
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

