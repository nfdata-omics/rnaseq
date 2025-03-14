#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("stringr"))
suppressMessages(library("openxlsx"))

### Script to merge results from gsea performed on different collections
option_list <- list(
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] list_of_result_tables",
                     option_list = option_list, prog = "gsea_merge",
                     description = "
                     Merge all the provided GSEA results into a single excel file.
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


excel_list = list()

for(f in files){
  # Filenames manipulation - The result may be a little dirty if the second group of the DGE comparison (e.g. 'B' in a comparison like Group_A_vs_B) contains a dot in its name.
  set_name = sub("^[^.]+\\.", "", strsplit(gsub(".xlsx", "", gsub("gsea.", "", basename(f))), "_vs_")[[1]][2])
  outname = gsub(".xlsx", "", gsub(set_name, "", basename(f)))
  
  # Importing result table
  dat = read.xlsx(f)
  # Adding the table as an additional sheet in the merged excel file
  excel_list[[set_name]] = dat
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