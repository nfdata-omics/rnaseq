#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("stringr"))
suppressMessages(library("openxlsx"))

### Script to merge results from gsea performed on different collections
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.1, help="The FDR cutoff to define significant pathways. [default \"%default\"]"),
  make_option(c("-n", "--num"), action="store", type="integer", default=25, help="The number of top most significantly enriched pathways to extract. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] list_of_result_tables",
                     option_list = option_list, prog = "gsea_merge",
                     description = "
                     Merge all the provided GSEA results into a single excel file, and produce a summary .txt file with the top N most enriched pathways.
                     list_of_result_tables must be the list of the complete paths to all the GSEA results to consider.
                     Multiple paths must be separated by blank-spaces, e.g. path/to/file1 path/to/file2 path/to/file3")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments=TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)<1 & !version) {
  stop("At least one result table must be supplied (the file containing GSEA results)", call.=FALSE)
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


files = arguments$args[1:length(arguments$args)]
fdr = opt$FDR
num = opt$num


excel_list = list()
sign_path = c()

for(f in files){
	if (grepl("_vs_", f)) {
		# Filenames manipulation for DGE comparisons - The result may be a little dirty if the second group of the DGE comparison (e.g. 'B' in a comparison like Group_A_vs_B) contains a dot in its name.
		set_name = sub("^[^.]+\\.", "", strsplit(gsub("\\.xlsx$", "", gsub("gsea.", "", basename(f))), "_vs_")[[1]][2])
		outname = gsub("\\.xlsx$", "", gsub(set_name, "", basename(f)))
	} else {
		# Filenames manipulation for DGE on continuous variables - The result may be a little dirty if the continuous variables contains a dot in its name.
		outname = paste0(paste(strsplit(gsub("\\.xlsx$", "", basename(f)), "\\.")[[1]][1:2], collapse="."), ".")
		set_name = gsub("\\.xlsx$", "", gsub(outname, "", basename(f)))
	}

  # Importing result table
  dat = read.xlsx(f)
  # Adding the table as an additional sheet in the merged excel file
  excel_list[[set_name]] = dat

  # Extracting and merging significant pathways
  if(!"empty"%in%colnames(dat)){
    dat$collection = set_name
    sign_path = rbind(sign_path, dat[dat$qvalue<fdr,])
  }
}


# Extracting the top N most significantly enriched pathways
if(length(sign_path)>0){
	sign_path = sign_path[order(sign_path$qvalue),]
	top_n = sign_path[1:min(num, nrow(sign_path)),]
	# Reshaping top_n in a convenient format for multiQC embedding
	top_n$minus.log.padj = -log10(top_n$qvalue)
	top_n = top_n[,c("ID","NES","pvalue","p.adjust","qvalue","minus.log.padj","collection")]
	colnames(top_n) = c("Pathway","NES","pvalue","p.adjust","qvalue","minus.log.padj","collection")
	write.table(top_n, paste(outname, "TOP_", num, ".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
} else {
	no_enrichment = data.frame(empty="No significant enrichment was found.")
	write.table(x=no_enrichment, file=paste(outname, "TOP_", num, ".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
}


# Saving the merged excel file, with different collections on different sheets
if(length(excel_list)>0){
  write.xlsx(excel_list, file=paste(outname, "xlsx", sep=""))
}




w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()
