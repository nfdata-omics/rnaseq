#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("sva"))
suppressMessages(library("stringr"))
suppressMessages(library("SummarizedExperiment"))

### Script to perform batch correction using the ComBat_seq function from package sva: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf 

option_list <- list(
  make_option(c("-k", "--keep"), action="store", type="character", default="", help="The name(s) of the metadata column(s) containing biological groups to preserve (signals from these variables are kept in data after adjustment). If multiple names are provided, they must be comma-separated with no blank spaces (e.g. variable1,variable2). Only categorical variables could be included. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] rna_object batch",
                     option_list = option_list, prog = "combat_batch_adj",
                     description = "Perform batch effect adjustment with the ComBat_seq algorithm implemented in the sva R package.
                     'rna_object' is the path of a .rds object containing a Summarized Experiment with expression data and samples metadata.
                     'batch' is the name of the metadata column containing batches. Only one batch variable could be adjuster for."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (rna_object and batch)", call.=FALSE)
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
batch = arguments$args[2]
biol_groups = str_split_1(opt$keep,",")



# Importing rna_object and extracting raw counts
rna_exp = get(load(rna_object))
exprs_to_adj = as.matrix(rna_exp@assays@data$counts)

# Metadata
meta = as.data.frame(colData(rna_exp))

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

# Updating the object
rna_exp@assays@data$combat_adjusted = exprs_to_adj
save(rna_exp, file=paste(substr(rna_object,1,nchar(rna_object)-4), ".batch_adj.rds", sep=""))

exprs_adj = cbind(rownames(exprs_adj), exprs_adj)
colnames(exprs_adj)[1] = "gene_name"

# Since the output is formatted as an adjusted count matrix with integer values, I re-attach gene information if available
if(rna_exp@metadata$gene_length){
  exprs_adj = cbind(exprs_adj[,1], rowData(rna_exp), exprs_adj[,2:ncol(exprs_adj)])
  colnames(exprs_adj)[1] = "gene_name"
}

write.table(exprs_adj, "counts_batch_adj.txt", quote=F, row.names=F, col.names=T, sep="\t")







########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


