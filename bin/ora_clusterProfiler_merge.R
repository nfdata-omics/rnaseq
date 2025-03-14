#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("stringr"))
suppressMessages(library("openxlsx"))

### Script to merge results from clusterProfiler over-representation analysis and to extract top results
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.1, help="The FDR cutoff to define significant pathways. [default \"%default\"]"),
  make_option(c("-n", "--num"), action="store", type="integer", default=25, help="The number of top most significantly enriched pathways to extract. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] list_of_result_tables",
                     option_list = option_list, prog = "ora_clusterProfiler_merge",
                     description = "
                     Produce an excel file containing all the significant pathways extracted from the provided list of clusterProfiler results.
                     It also produces a reduced .txt file containing the top N most significantly enriched pathways.
                     list_of_result_tables must be the list of the complete paths to all the ORA-clusterProfiler results to consider.
                     Multiple paths must be separated by blank-spaces, e.g. path/to/file1 path/to/file2 path/to/file3")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments=TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)<1 & !version) {
  stop("At least one result table must be supplied (the files containing clusterProfiler results)", call.=FALSE)
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



# Filenames manipulation
split = strsplit(files[1], "\\.")[[1]]
suffix = paste("", split[length(split)-1], split[length(split)], sep=".")

all_sign = c()
excel_list = list()

for(f in files){
  # Filenames manipulation - The result may be a little dirty if the second group of the DGE comparison (e.g. 'B' in a comparison like Group_A_vs_B) contains a dot in its name.
  set_name = sub("^[^.]+\\.", "", strsplit(gsub(suffix, "", gsub("ora_CP.", "", basename(f))), "_vs_")[[1]][2])
  outname = gsub(suffix, "", gsub(set_name, "", basename(f)))
  
  # Importing result table
  dat = read.delim(f, h=T)
  
  # If the table is not empty, extract the significant records
  if(colnames(dat)[1] != "empty"){
    signif = dat[dat$p.adjust < fdr,]
    if(dim(signif)[1] > 0){
      excel_list[[set_name]] = signif
      signif$collection = set_name
      all_sign = rbind(all_sign, signif)
    }
  }
  
}


# Saving top N most enriched pathways, among the significant ones
if(!is.null(nrow(all_sign))){
  all_sign = all_sign[order(all_sign$p.adjust),]
  top_n = all_sign[1:min(num, nrow(all_sign)),]
  # Reshaping top_n in a convenient format for multiQC embedding
  top_n$minus.log.padj = -log10(top_n$p.adjust)
  top_n = top_n[,c(1,5,8,9,14,13)]
  colnames(top_n) = c("Pathway", "gene.ratio", "p.value", "p.adjust", "minus.log.padj", "collection")
  write.table(top_n, paste(outname, "TOP_", num, suffix, sep=""), col.names=T, row.names=F, quote=F, sep="\t")
}


# Saving all significant result to an excel table, with different collections on different sheets
if(length(excel_list)>0){
  write.xlsx(excel_list, file=paste(outname, substr(suffix, 2, nchar(suffix)-3), "xlsx", sep=""))
}





########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


