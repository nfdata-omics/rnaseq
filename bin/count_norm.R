#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("edgeR"))
suppressMessages(library("limma"))
suppressMessages(library("SummarizedExperiment"))

### Script to perform raw counts normalization into cpm and/or rpkm

option_list <- list(
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] rna_object",
                     option_list = option_list, prog = "count_norm",
                     description = "Perform raw counts normalization into cpm and possibly rpkm. Rpkm are calculated only if raw counts were originally produced using featureCounts.
                     'rna_object' is the path of a .rds object containing a Summarized Experiment with raw counts and samples metadata."
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


# Importing rna_object
rna_exp = get(load(rna_object))

# Creating edgeR object
y_all = DGEList(counts=rna_exp@assays@data$counts, genes=rowData(rna_exp))

# Cpm
cpm = cpm(y_all, log=FALSE)
cpm = cbind(rownames(cpm), cpm)
colnames(cpm)[1] = "gene_name"
write.table(cpm, "cpm.txt", quote=F, row.names=F, col.names=T, sep="\t")
rna_exp@assays@data$cpm = cpm(y_all, log=FALSE)

if(rna_exp@metadata$featureCounts){
	
	# Log-norm-rpkm (TMM normalization)
	y_all = calcNormFactors(y_all, method="TMM")
	log_norm_rpkm = rpkm(y_all, log=T, gene.length=y_all$genes$Length)
	log_norm_rpkm = cbind(rownames(log_norm_rpkm), log_norm_rpkm)
	colnames(log_norm_rpkm)[1] = "gene_name"
	write.table(log_norm_rpkm, "log.norm.rpkm.txt", quote=F, row.names=F, col.names=T, sep="\t")
	rna_exp@assays@data$log_norm_rpkm = rpkm(y_all, log=T, gene.length=y_all$genes$Length)
	
	# Lib size and normalization factor
	size_factors = data.frame(samples=rownames(y_all$samples), lib.size=y_all$samples$lib.size, norm.factors=y_all$samples$norm.factors)
	write.table(size_factors, "lib_size_factors.txt", quote=F, row.names=F, col.names=T, sep="\t")

} else {

	# Lib size
	size_factors = data.frame(samples=rownames(y_all$samples), lib.size=y_all$samples$lib.size)
	write.table(size_factors, "lib_size_factors.txt", quote=F, row.names=F, col.names=T, sep="\t")

}

# Updating the object
save(rna_exp, file=paste(substr(rna_object,1,nchar(rna_object)-4),".norm.rds",sep=""))




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
 print(w)
}
sink()



