#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))

### Script to fit a DGE model using DESeq2
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.05, help="The FDR cutoff to define significant genes. [default \"%default\"]"),
  make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output filenames, e.g. deseq2_toptable.variable_num_vs_denom_SUFFIX.txt. [default \"%default\"]")
)
### change logFC names
parser<-OptionParser(usage = "%prog [options] input_model contrast",
                     option_list = option_list, prog = "dge_deseq2_results",
                     description = "Extract results from a DGE model fitted with DESeq2. 
                     'input_model' is the path of the DESeq2 object containing the already fitted model.
                     'contrast' is a string of format variable:test:reference, with words separated by colon (:),
                                where variable is the name of the variable (metadata column) on which the comparison will be performed (e.g. genotype - It must be one of the variables previously included in the model formula.)
                                test is the name of the factor level to be tested (the numerator of the log2 Fold Change, e.g. KO)
                                reference is the name of the factor level to be used as reference (the denominator of the log2 Fold Change, e.g. WT)."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=2) {
  stop("Two arguments must be supplied (input_model and contrast)", call.=FALSE)
}

input_model = arguments$args[1]
contrast = str_split_1(arguments$args[2], ":")
variable = contrast[1]
num = contrast[2]
denom = contrast[3]
fdr = opt$FDR
suffix = opt$suffix

# Importing dge object
dds = get(load(input_model))

# Extracting results
my_res = results(dds, contrast=c(variable, num, denom), independentFiltering=TRUE, cooksCutoff=FALSE, alpha=fdr)
# I turn off cooksCutoff for outlier detection, but in the model fitting function there was minRepforReplace=7
my_res = my_res[order(my_res$pvalue, decreasing=F),]
capture.output(summary(my_res, alpha=fdr), file=paste("dge_summary.",variable,"_",num,"_vs_",denom,suffix,".txt", sep=""))

pdf(paste("MAplot.",variable,"_",num,"_vs_",denom,".pdf", sep=""))
DESeq2::plotMA(my_res, ylim=c(-6,6), alpha=fdr, main=paste(num," vs ",denom, sep=""))
dev.off()

# Saving toptable
toptable = my_res[,c("baseMean","log2FoldChange","pvalue","padj")]
toptable = cbind(rownames(toptable), toptable)
colnames(toptable)[1] = "gene_name"
write.table(toptable, paste("deseq2_toptable.",variable,"_",num,"_vs_",denom,suffix,".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")


# Saving package versions
x = sessionInfo()
my_pkgs = c()
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

