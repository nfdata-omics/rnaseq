#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("edgeR"))
suppressMessages(library("limma"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))
suppressMessages(library("SummarizedExperiment"))

### Script to fit a DGE model using DESeq2
option_list <- list(
  #make_option(c("-n", "--num_samples"), action="store", type="integer", default=1, help="Number of samples to define expressed genes. Expressed genes are defined as those genes showing at least 1 cpm on at least 'n' different samples. Usually, this is the size of the smallest sample group. [default \"%default\"]"),
  make_option(c("-e", "--expressed"), action="store", type="double", default=0.5, help="Fraction of samples required to define a gene as expressed. A gene is considered expressed if it shows at least 1 cpm in at least a fraction 'e' of the samples. [default \"%default\"]"),
  make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output paths, e.g. deseq2_obj_SUFFIX.rds. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] rna_object model_formula",
                     option_list = option_list, prog = "dge_deseq2_fit",
                     description = "Fit a differential gene expression (DGE) model using DESeq2. 
                     'rna_object' is the path of a .rds object containing a Summarized Experiment with expression data and samples metadata.
                     'model_formula' is the formula used for the design of the DESeq2 model. It must start with a '~' and include all the biological and technical variables that should be accounted for (e.g. ~genotype+treatment+batch) with no blank spaces."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (rna_object and model_formula)", call.=FALSE)
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
model_formula = as.formula(arguments$args[2])
expr = opt$expressed
suffix = opt$suffix


# Importing rna_object
rna_exp = get(load(rna_object))
counts = rna_exp@assays@data$counts

# Defining expressed genes
num_samples = round(expr*nrow(colData(rna_exp)))
isexpr = rowSums(rna_exp@assays@data$cpm>1) >= num_samples

# Defining the DESeq2 object
dds = DESeqDataSetFromMatrix(countData=as.matrix(counts[isexpr,]), colData=colData(rna_exp), design=model_formula)

# DispEst plot for the whole dataset
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
pdf(paste("DispEst_plot",suffix,".pdf", sep=""), width=8, height=8)
plotDispEsts(dds)
dev.off()

# Fitting models
dds = DESeq(object=dds, test="Wald", fitType="parametric", betaPrior=FALSE, minReplicatesForReplace=7, parallel=F)
save(dds, file=paste("deseq2_obj",suffix,".rds", sep=""))



########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


