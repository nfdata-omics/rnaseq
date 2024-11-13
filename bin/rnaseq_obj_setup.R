#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("edgeR"))
suppressMessages(library("limma"))
suppressMessages(library("SummarizedExperiment"))

### Script to setup the summarized experiment object containing raw counts and metadata

option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix has been produced using featureCounts or not. If true, gene names are supposed to be reported in the first column, and gene informations (Chr/Start/End/Strand/Length) are supposed to be reported in the 2-6 columns. If false, only the first column containing gene names is expected. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

################################
# Consider to add an option/parameter to specify the output name! In this way batch-corrected data (for instance) could be saved in a different object

parser<-OptionParser(usage = "%prog [options] raw_counts metadata",
                     option_list = option_list, prog = "rnaseq_obj_setup",
                     description = "Setup the summarized experiment object containing raw counts and sample metadata, that will be used as input for the subsequent pipeline steps. 
                     'raw_counts' is the path of a tab delimited file containing the raw counts matrix.
                     'metadata' is the path of a csv file containing sample metadata. Sample names must be reported in the first column. Columns containing numbers are assumed to be quantitative variables."
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
feature_counts = opt$featureCounts

# Importing raw counts
counts = read.delim(raw_counts, h=T, row.names=1, check.names=F)

# Importing metadata
meta = read.csv(metadata, h=T, row.names=1, check.names=F)
# Numeric variables are assumed to be quantitative, other variables are assumed to be categorical
for(i in 1:ncol(meta)){
  if(!is.numeric(meta[,i])){
    meta[,i] = as.factor(meta[,i])
  }
}

# Creating the summarized exp object
if(feature_counts){
  
  # Matching names
  meta = meta[rownames(meta)%in%colnames(counts)[6:ncol(counts)],]
  meta = meta[match(colnames(counts)[6:ncol(counts)], rownames(meta)),]
  # Creating the summarized experiment object
  rna_exp = SummarizedExperiment(assays=list(counts=as.matrix(counts[,6:ncol(counts)])), rowData=as.data.frame(counts[,1:5]), colData=meta)
  metadata(rna_exp)$featureCounts = TRUE
  
} else {
  
  # Matching names
  meta = meta[rownames(meta)%in%colnames(counts),]
  meta = meta[match(colnames(counts), rownames(meta)),]
  # Summ exp obj
  rna_exp = SummarizedExperiment(assays=list(counts=as.matrix(counts)), rowData=as.data.frame(rownames(counts)), colData=meta)
  metadata(rna_exp)$featureCounts = FALSE
  
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


