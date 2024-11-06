#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("edgeR"))
suppressMessages(library("limma"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))

### Script to fit a DGE model using DESeq2
option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix has been produced using featureCounts or not. If true, gene names are supposed to be reported in the first column, and gene informations (Chr/Start/End/Strand/Length) are supposed to be reported in the 2-6 columns. [default \"%default\"]"),
  make_option(c("-n", "--num_samples"), action="store", type="integer", default=1, help="Number of samples to define expressed genes. Expressed genes are defined as those genes showing at least 1 cpm on at least 'n' different samples. Usually, this is the size of the smallest sample group. [default \"%default\"]"),
  make_option(c("-c", "--categorical"), action="store", type="character", default="", help="The name(s) of categorical variables included in the model. They must be the names of the corresponding metadata columns. If multiple names are provided, they must be comma-separated with no blank spaces (e.g. genotype,treatment). The other numerical variables in the model_formula are assumed to be quantitative. [default \"%default\"]"),
  make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output paths, e.g. deseq2_obj_SUFFIX.rds. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] input_data metadata model_formula",
                     option_list = option_list, prog = "dge_deseq2_fit",
                     description = "Fit a differential gene expression (DGE) model using DESeq2. 
                     'input_data' is the path of a tab delimited file containing the raw counts matrix.
                     'metadata' is the path of a tab delimited file containing sample metadata. Sample names must be reported in the first column.
                     'model_formula' is the formula used for the design of the DESeq2 model. It must start with a '~' and include all the biological and technical variables that should be accounted for (e.g. ~genotype+treatment+batch) with no blank spaces."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=3) {
  stop("Three arguments must be supplied (input_data)", call.=FALSE)
}

input_data = arguments$args[1]
metadata = arguments$args[2]
model_formula = as.formula(arguments$args[3])
feature_counts = opt$featureCounts
categorical_vars = str_split_1(opt$categorical, ",")
num_samples = opt$num_samples
suffix = opt$suffix


# Importing raw counts
counts = read.delim(input_data, h=T, row.names=1, check.names=F)
if(feature_counts){
  counts = counts[,6:ncol(counts)]
}
# Defining expressed genes
y_all <- DGEList(counts=counts)
isexpr = rowSums(cpm(y_all)>1) >= num_samples


# Importing metadata and matching names
meta = read.delim(metadata, h=T, row.names=1, check.names=F)
meta = meta[rownames(meta)%in%colnames(counts),]
meta = meta[match(colnames(counts), rownames(meta)),]
# Converting categorical variables into factors
if(categorical_vars[1]!="" & sum(!categorical_vars%in%colnames(meta))>0){
  stop("Categorical variables must be names of metadata columns. If more than one variable are provided, they must be comma-separated with no blank spaces (e.g. treatment,group)")
}
for(i in 1:ncol(meta)){
  if(colnames(meta)[i]%in%categorical_vars)
    meta[,i] = as.factor(meta[,i])
}


# Defining the DESeq2 object
dds = DESeqDataSetFromMatrix(countData=as.matrix(counts[isexpr,]), colData=meta, design=model_formula)

# DispEst plot for the whole dataset
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
pdf(paste("DispEst_plot",suffix,".pdf", sep=""), width=8, height=8)
plotDispEsts(dds)
dev.off()

# Fitting models
dds = DESeq(object=dds, test="Wald", fitType="parametric", betaPrior=FALSE, minReplicatesForReplace=7, parallel=F)
save(dds, file=paste("deseq2_obj",suffix,".rds", sep=""))




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


