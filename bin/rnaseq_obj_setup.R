#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("SummarizedExperiment"))

### Script to setup the summarized experiment object containing raw counts and metadata

option_list <- list(
  make_option(c("-c","--col_num"), action="store", type="integer", default=1, help="Number of columns in the raw_counts matrix containing gene information (e.g. GeneID, Chr, Start, End, Strand, Length). These columns must be the first ones of the matrix. If gene lengths are provided, the corresponding column must be called 'Length'. [default \"%default\"]"),
  make_option(c("-i","--gene_id"), action="store", type="integer", default=1, help="The column in the raw_counts matrix that contains the gene name identifiers. Gene identifiers should be unique; if they are not, a sequential number will be automatically added to make each identifier unique. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

################################
# Consider to add an option/parameter to specify the output name! In this way batch-corrected data (for instance) could be saved in a different object

parser<-OptionParser(usage = "%prog [options] raw_counts metadata",
                     option_list = option_list, prog = "rnaseq_obj_setup",
                     description = "Setup the summarized experiment object containing raw counts and sample metadata, that will be used as input for the subsequent pipeline steps. 
                     'raw_counts' is the path of a tab delimited file containing the raw counts matrix.
                     'metadata' is the path of a csv file containing sample metadata. Sample names must be reported in the first column. Columns containing numbers are assumed to be quantitative variables.
                     .META: raw_counts
                        1 ... c : gene information (e.g. GeneID, position, strand, Length, ...)
                        i (<= c): gene IDs
                        c+1 ... : expression values"
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (raw_counts and metadata)", call.=FALSE)
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

raw_counts = arguments$args[1]
metadata = arguments$args[2]
col_num = opt$col_num
id = opt$gene_id

if(id > col_num){
  stop("The column containing gene identifiers must be one of the first #col_num columns!")
}

# Importing raw counts
counts = read.delim(raw_counts, h=T)
# If gene identifiers are not unique, they will be forced to be
if(length(unique(counts[,id]))<length(counts[,id])){
  counts[,id] = make.unique(counts[,id])
}
rownames(counts) = counts[,id]
counts = counts[,-id]

# Importing metadata
meta = read.csv(metadata, h=T, row.names=1, check.names=F)
# Numeric variables are assumed to be quantitative, other variables are assumed to be categorical
for(i in 1:ncol(meta)){
  if(!is.numeric(meta[,i])){
    meta[,i] = as.factor(meta[,i])
  }
}

# Matching names
meta = meta[rownames(meta)%in%colnames(counts)[(col_num):ncol(counts)],]
meta = meta[match(colnames(counts)[(col_num):ncol(counts)], rownames(meta)),]

# Creating the summarized exp object
if(col_num==1){
  
  rna_exp = SummarizedExperiment(assays=list(counts=as.matrix(counts)), rowData=as.data.frame(rownames(counts)), colData=meta)
  colnames(rowData(rna_exp)) = "gene_name"
  metadata(rna_exp)$gene_length = FALSE
  
} else {
  
  rna_exp = SummarizedExperiment(assays=list(counts=as.matrix(counts[,(col_num):ncol(counts)])), rowData=as.data.frame(counts[,1:(col_num-1)]), colData=meta)
  colnames(rowData(rna_exp)) = colnames(counts)[1:(col_num-1)]
  # Is gene length information available?
  if("Length"%in%colnames(counts)[1:(col_num-1)]){
    metadata(rna_exp)$gene_length = TRUE
  } else {
    metadata(rna_exp)$gene_length = FALSE
  }
  
}



# Saving the object
save(rna_exp, file="rna_SummExp.rds")




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


