#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("pheatmap"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("SummarizedExperiment"))

### Script to extract expression values of specific control genes 

option_list <- list(
  make_option(c("-a", "--annotations"), action="store", type="character", default="", help="The name(s) of categorical variables to highlight on the top of the heatmap (as annotation column bars). They must be the names of the corresponding metadata columns. If multiple names are provided, they must be comma-separated with no blank spaces (e.g. genotype,treatment). [default \"%default\"]"),
  make_option(c("-d", "--dendrogram"), action="store", type="integer", default=2, help="Whether to perform clustering and show the corresponding dendrogram on rows/genes ('1'), on columns/samples ('2'), both ('3') or none ('0'). [default \"%default\"]"),
  make_option(c("-s", "--scale_rows"), action="store_true", default=FALSE, help="Whether the values should be centered and scaled in the row direction. If true, genes with constant expression will be removed. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] rna_object genes_list
                     'rna_object' is the path of a .rds object containing a Summarized Experiment with expression data and samples metadata.
                     'genes_list' is the path to a tab-delimited file containing the list of control genes to evaluate. Each gene must be reported in a different row, with no header. Gene encoding must be the same used for the raw counts matrix.",
                     option_list = option_list, prog = "control_genes",
                     description = "Extract the expression values (cpm or rpkm, depending on whether gene length information was originally available) of specific control genes, and produce the corresponding heatmap."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (rna_object and genes_list)", call.=FALSE)
}

if(version){
  # Printing package versions
  x = sessionInfo()
  cat(paste(" "," "," R: ", x$R.version$major,".", x$R.version$minor, "\n", sep=""))
  for(i in 1:length(x$otherPkgs)){
    cat(paste(" "," "," ", x$otherPkgs[[i]]$Package,": ", x$otherPkgs[[i]]$Version, "\n", sep=""))
  }
  quit()
}

rna_object = arguments$args[1]
genes_list = arguments$args[2]
annot = str_split_1(opt$annotations, ",")
dend = opt$dendrogram
scale_rows = opt$scale_rows


# Importing rna_object
rna_exp = get(load(rna_object))

# Extracting expression data - cpm or rpkm, depending on whether gene length information was available
if(rna_exp@metadata$gene_length){
  data = rna_exp@assays@data$log_norm_rpkm
} else {
  data = rna_exp@assays@data$cpm
}

# Importing control genes
control_genes = as.character(read.delim(genes_list, h=F)$V1)
control_genes = control_genes[control_genes%in%rownames(data)]
if(length(control_genes)<=1){
  stop("Only one or zero control genes were retrieved in the expression matrix. Check consistency between the expression data and the supplied list of genes.")
}

# Filtering data table
data_ctrl = data[control_genes,]

# Excluding constant genes if scaling on rows
if(scale_rows){
  var = apply(data_ctrl, 1, var)
  data_ctrl = data_ctrl[var!=0,]
  if(is.null(dim(data_ctrl))){
    stop("All the selected genes have a constant expression across the samples and were excluded. 
         (It is possible that one single gene has a non-constant expression, but heatmap for a single gene is not supported.)")}
  scale_value = "row"
} else {
  scale_value = "none"
}

# Annotations for heatmap
meta = colData(rna_exp)
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

png("heatmap_control_genes_mqc.png", width=8, height=10, units="in", res=1200)
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

