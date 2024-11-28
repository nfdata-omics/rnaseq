#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("limma"))
suppressMessages(library("edgeR"))
suppressMessages(library("SummarizedExperiment"))

### Script to perform PCA and MDS

option_list <- list(
  #make_option(c("-n", "--num_samples"), action="store", type="integer", default=1, help="Number of samples to define expressed genes. Expressed genes are defined as those genes showing at least 1 cpm on at least 'n' different samples. Usually, this is the size of the smallest sample group. [default \"%default\"]"),
  make_option(c("-e", "--expressed"), action="store", type="double", default=0.5, help="Fraction of samples required to define a gene as expressed. A gene is considered expressed if it shows at least 1 cpm in at least a fraction 'e' of the samples. [default \"%default\"]"),
  make_option(c("-t", "--top_var"), action="store", type="integer", default=5000, help="Number of the most variable genes to consider. PCA and MDS will be performed on the top 't' most variable genes, in order to focus on the main sources of variability in the dataset. [default \"%default\"]"),
  make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output paths, e.g. PCA_scores_SUFFIX.txt. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] rna_object",
                     option_list = option_list, prog = "pca_mds",
                     description = "Perform PCA and MDS on the expression values (cpm or rpkm for PCA and limma-voom normalized counts for MDS) of the most variable genes among the expressed ones, as defined by the options.
                     'rna_object' is the path of a .rds object containing a Summarized Experiment with expression data and samples metadata."
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
expr = opt$expressed
top_var = opt$top_var
suffix = opt$suffix


# Importing rna_object
rna_exp = get(load(rna_object))

# Importing expression values
if(rna_exp@metadata$gene_length){
  data = rna_exp@assays@data$log_norm_rpkm
  cpm = rna_exp@assays@data$cpm
} else {
  data = rna_exp@assays@data$cpm
  cpm = data
}

# Number of expressed genes
num_samples = round(expr*nrow(colData(rna_exp)))
isexpr = rowSums(cpm>1) >= num_samples
# Filtering not-expressed genes
data = data[isexpr,]

# Considering the top t most variable genes
var = apply(data, 1, var)
names(var) = rownames(data)
sorted_var = sort(var, decreasing=T)
data = data[names(sorted_var[1:top_var]),]
# Saving data to produce dendrogram directly from python
data2dendr = cbind(rownames(data), data)
colnames(data2dendr)[1] = "gene_name"
write.table(data2dendr, paste("exprs_data_to_hc_dendrog",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")


# Metadata
meta = colData(rna_exp)

# PCA
pca = prcomp(t(data), scale=T, center=T)
score = as.data.frame(pca$x)
score = cbind(rownames(score), score)
colnames(score)[1]="samples"
score = cbind(score, meta)
write.table(score, paste("PCA_scores",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")

# Saving variance for screeplot
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
var_df = data.frame(PC=paste("PC",seq(1:length(var)), sep=""), var=var)
write.table(var_df, paste("PCA_explained_variance",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")



# MDS
y_all = DGEList(counts=rna_exp@assays@data$counts, genes=rowData(rna_exp))
y_all = y_all[names(sorted_var[1:top_var]),]
voom = voom(y_all, plot=F)
mds = plotMDS(voom, main="MDS plot", plot=F)
mds_score = data.frame(samples=rownames(mds@.Data[[5]]), x=mds$x, y=mds$y)
meta = meta[match(mds_score$samples, rownames(meta)),] # probably useless but one more it's better than one less
mds_score = cbind(mds_score, meta)
write.table(mds_score, paste("MDS_scores",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")





########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()



