#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("sva"))
suppressMessages(library("stringr"))

### Script to perform batch correction using the ComBat_seq function from package sva: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf 

option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix was originally produced using featureCounts or not. [default \"%default\"]"),
  make_option(c("-k", "--keep"), action="store", type="character", default="", help="The name(s) of the metadata column(s) containing biological groups to preserve (signals from these variables are kept in data after adjustment). If multiple names are provided, they must be comma-separated with no blank spaces (e.g. variable1,variable2). Only categorical variables could be included. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] input_data metadata batch",
                     option_list = option_list, prog = "combat_batch_adj",
                     description = "Perform batch effect adjustment with the ComBat_seq algorithm implemented in the sva R package.
                     'input_data' is the path of a tab delimited file containing the raw counts matrix.
                     'metadata' is the path of a tab delimited file containing sample metadata. Sample names must be reported in the first column.
                     'batch' is the name of the metadata column containing batches. Only one batch variable could be adjuster for."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=3) {
  stop("Three arguments must be supplied (input_data, metadata and batch)", call.=FALSE)
}

input_data = arguments$args[1]
metadata = arguments$args[2]
batch = arguments$args[3]
feature_counts = opt$featureCounts
biol_groups = str_split_1(opt$keep,",")


# Importing raw counts
counts = read.delim(input_data, h=T, row.names=1, check.names=F)
# Considering only expression values (excluding gene informations when available)
if(feature_counts){
  exprs_to_adj = counts[,6:ncol(counts)]
  exprs_to_adj = as.matrix(exprs_to_adj)
} else {
  exprs_to_adj = counts
  exprs_to_adj = as.matrix(exprs_to_adj)
}

# Importing metadata and matching names
meta = read.delim(metadata, h=T, row.names=1, check.names=F)
meta = meta[rownames(meta)%in%colnames(exprs_to_adj),]
meta = meta[match(colnames(exprs_to_adj), rownames(meta)),]

# Handling batch variable
if(!batch%in%colnames(meta)){
  stop("The specified batch must be a column in the metadata table. Check the correct name.")
}
meta[,batch] = as.factor(meta[,batch])
to_adj = meta[,batch]

# Handling biological covariates
for(i in 1:length(biol_groups)){
  if(biol_groups[i]!="" & !biol_groups[i]%in%colnames(meta)){
    stop("The specified biological group(s) must be column(s) in the metadata table. Check the correct names. If more than one variable is specified, their names must be comma-separated withoud blank spaces.")
  }
}
if(biol_groups[1]!=""){
  to_keep = meta[,biol_groups] 
}

# Performing adjustment with ComBat_seq
if(biol_groups[1]==""){
  exprs_adj = ComBat_seq(exprs_to_adj, batch=to_adj, group=NULL)
}
if(biol_groups[1]!="" & length(biol_groups)==1){
  exprs_adj = ComBat_seq(exprs_to_adj, batch=to_adj, group=to_keep)
}
if(length(biol_groups)>1){
  exprs_adj = ComBat_seq(exprs_to_adj, batch=to_adj, covar_mod=to_keep)
}
exprs_adj = cbind(rownames(exprs_adj), exprs_adj)
colnames(exprs_adj)[1] = "gene_name"

# Since the output is formatted as an adjusted count matrix with integer values, I re-attach gene information if available
if(feature_counts){
  exprs_adj = cbind(exprs_adj[,1], counts[,1:5], exprs_adj[,2:ncol(exprs_adj)])
  colnames(exprs_adj)[1] = "gene_name"
}

write.table(exprs_adj, "counts_batch_adj.txt", quote=F, row.names=F, col.names=T, sep="\t")



# Saving package versions
x = sessionInfo()
my_pkgs = c(paste(" "," "," ","R: ",x$R.version$major,".",x$R.version$minor, sep=""))
for(i in 1:length(x$otherPkgs)){
  my_pkgs = c(my_pkgs, paste(" "," "," ", x$otherPkgs[[i]]$Package,": ", x$otherPkgs[[i]]$Version, sep=""))
}
pkgVersion = file("R_pkgs_versions.txt")
writeLines(my_pkgs, pkgVersion)
close(pkgVersion)




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


