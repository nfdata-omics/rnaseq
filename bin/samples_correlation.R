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
  make_option(c("-e", "--expressed"), action="store", type="double", default=0.5, help="Fraction of samples required to define a gene as expressed. A gene is considered expressed if it shows at least 1 cpm in at least a fraction 'e' of the samples. [default \"%default\"]"),
  make_option(c("-s", "--spearman"), action="store_true", default=FALSE, help="Whether Spearman's correlation should be calculated instead of Pearsons' correlation. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] rna_object
                     'rna_object' is the path of a .rds object containing a Summarized Experiment with expression data and samples metadata.",
                     option_list = option_list, prog = "samples_correlation",
                     description = "Calculate the correlation between samples based on their normalized expression values, and produce the corresponding heatmap."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=1 & !version) {
  stop("One argument must be supplied (rna_object)", call.=FALSE)
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
annot = str_split_1(opt$annotations, ",")
expr = opt$expressed
spearman = opt$spearman


# Importing rna_object
rna_exp = get(load(rna_object))
cpm = rna_exp@assays@data$cpm

# Number of expressed genes
num_samples = round(expr*nrow(colData(rna_exp)))
isexpr = rowSums(cpm>1) >= num_samples
# Filtering not-expressed genes
cpm = cpm[isexpr,]


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


# Correlation
if(spearman){
  my_cor = cor(cpm, method="spearman")
  my_cor_table = cbind(rownames(my_cor), my_cor)
  colnames(my_cor_table)[1] = "Spearman_cor"
} else {
  my_cor = cor(cpm)
  my_cor_table = cbind(rownames(my_cor), my_cor)
  colnames(my_cor_table)[1] = "Pearson_cor"
}
write.table(my_cor_table, "samples_correlation_table.txt", quote=F, row.names=F, col.names=T, sep="\t")


# Heatmap
pdf("samples_correlation_heatmap.pdf", height=10, width=12)
if(nrow(meta)<=20){
  pheatmap(my_cor, scale="none", display_numbers=T, treeheight_row=0, annotation_col=col_annot)
} else {
  pheatmap(my_cor, scale="none", display_numbers=F, treeheight_row=0, annotation_col=col_annot)
}
dev.off()







########################

# w=warnings()
# sink(stderr())
# if(!is.null(w)){
#   print(w)
# }
# sink()