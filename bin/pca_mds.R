#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("limma"))
suppressMessages(library("edgeR"))

### Script to perform PCA and MDS

option_list <- list(
  make_option(c("-f", "--featureCounts"), action="store_true", default=FALSE, help="Whether the raw counts matrix was originally produced using featureCounts or not. If true, PCA will be performed on rpkm values instead of cpm. [default \"%default\"]"),
  make_option(c("-n", "--num_samples"), action="store", type="integer", default=1, help="Number of samples to define expressed genes. Expressed genes will be defined as those genes showing at least 1 cpm on at least 'n' different samples. Usually, this is the size of the smallest sample group. [default \"%default\"]"),
  make_option(c("-t", "--top_var"), action="store", type="integer", default=5000, help="Number of the most variable genes to consider. PCA and MDS will be performed on the top 't' most variable genes, in order to focus on the main sources of variability in the dataset. [default \"%default\"]"),
 make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output paths, e.g. PCA_scores_SUFFIX.txt. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] metadata",
                     option_list = option_list, prog = "pca_mds",
                     description = "Perform PCA and MDS on the expression values (cpm or rpkm and limma-voom) of the most variable genes among the expressed ones, as defined by the options.
                     'metadata' is the path of a tab delimited file containing sample metadata. Sample names must be reported in the first column."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=1) {
  stop("One argument must be supplied (sample metadata)", call.=FALSE)
}

metadata = arguments$args[1]
feature_counts = opt$featureCounts
min_size = opt$num_samples
top_var = opt$top_var
suffix = opt$suffix

# Importing expression values
if(feature_counts){
  data = read.delim("log.norm.rpkm.txt", h=T, row.names=1)
  cpm = read.delim("cpm.txt", h=T, row.names=1)
} else {
  data = read.delim("cpm.txt", h=T, row.names=1)
  cpm = data
}

# Number of expressed genes
isexpr = rowSums(cpm>1) >= min_size
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
#write.table(data2dendr, "exprs_data_to_hc_dendrog.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(data2dendr, paste("exprs_data_to_hc_dendrog",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")


# Importing metadata and matching names
meta = read.delim(metadata, h=T, row.names=1, check.names=F)
meta = meta[rownames(meta)%in%colnames(data),]
meta = meta[match(colnames(data), rownames(meta)),]


# PCA
pca = prcomp(t(data), scale=T, center=T)
score = as.data.frame(pca$x)
score = cbind(rownames(score), score)
colnames(score)[1]="samples"
score = cbind(score, meta)
#write.table(score, "PCA_scores.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(score, paste("PCA_scores",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")

# Saving variance for screeplot
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
var_df = data.frame(PC=paste("PC",seq(1:length(var)), sep=""), var=var)
write.table(var_df, paste("PCA_explained_variance",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")



# MDS
y_all = get(load("dge_obj.rds"))
y_all = y_all[names(sorted_var[1:top_var]),]
voom = voom(y_all, plot=F)
mds = plotMDS(voom, main="MDS plot", plot=F)
mds_score = data.frame(samples=rownames(mds@.Data[[5]]), x=mds$x, y=mds$y)
meta = meta[match(mds_score$samples, rownames(meta)),] # probably useless but one more it's better than one less
mds_score = cbind(mds_score, meta)
#write.table(mds_score, "MDS_scores.txt", quote=F, row.names=F, col.names=T, sep="\t")
write.table(mds_score, paste("MDS_scores",suffix,".txt",sep=""), quote=F, row.names=F, col.names=T, sep="\t")




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



