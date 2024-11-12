#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))

### Script to extract the most significantly enriched pathways from all the enrichR results.
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.05, help="The FDR cutoff to define significant pathways. [default \"%default\"]"),
  make_option(c("-n", "--num"), action="store", type="integer", default=25, help="The number of top most significantly enriched pathways to extract. [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] enrichr_full_output",
                     option_list = option_list, prog = "ora_top_results",
                     description = "Extract the top N most significantly enriched pathways from all the collections previously tested with enrichR. 
                     'enrichr_full_output' is the path of the .xlsx file containing all the enrichR results from which significant pthways are extracted.")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=1) {
  stop("One argument must be supplied (enrichr_full_output)", call.=FALSE)
}

file = arguments$args[1]
fdr = opt$FDR
num = opt$num


# Alternatively, it is possible to not specify the input file (enrichr_full_output) and recursively apply the script to all the enrichr.*.xlsx tables in the current directory
# Listing all the output .xlsx tables previously produced with enrichR
#lf = list.files(pattern=glob2rx("enrichr.*.xlsx"))

# Defining a function that calculates the gene ratio from the Overlap statistic (e.g overlap=53/100 --> gene.ratio=0.53)
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))

#for(file in lf){
  # Creating a single table for each enrichment analysis (merging all the pathway collections)
  s = openxlsx::getSheetNames(file)
  enrichR.table = data.frame()
  for(dat in s){
    Table <- openxlsx::read.xlsx(xlsxFile=file, sheet=dat, startRow=1, colNames=T, rowNames=T,
                                 detectDates=F, skipEmptyRows=T, skipEmptyCols=T, na.strings="NA", fillMergedCells=F)
    Table$Collection = rep(dat, nrow(Table))
    enrichR.table = rbind(enrichR.table, Table)
  }
  
  # Considering only significant pathways
  p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < fdr,])
  
  if(length(p)>0) { # Discarding the not-significant results (to avoid errors)
    pathways.dataframe = data.frame(Pathway=p, gene.ratio=sapply(p, fx), p.value=enrichR.table[p,]$P.value, p.value.adj=enrichR.table[p,]$Adjusted.P.value, 
                                    combined.score=enrichR.table[p,]$Combined.Score, collection=enrichR.table[p,]$Collection)
    # Sorting and saving the dataframe with top results
    pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
    pathways.dataframe = pathways.dataframe[1:min(num,nrow(pathways.dataframe)),]
    write.table(pathways.dataframe, paste(str_sub(file,1,-6),'_TOP',num,'.txt',sep=''), row.names=F, col.names=T, quote=F, sep="\t")
  }
#}


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



