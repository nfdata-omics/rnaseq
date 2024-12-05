#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("enrichplot"))
#suppressMessages(library("msigdbr"))
suppressMessages(library("openxlsx"))

### Script to perform first line functional enrichment using GSEA from clusterProfiler
option_list <- list(
  make_option(c("-r", "--ranking"), action="store", type="character", default="logFC", help="The ranking statistics used to order genes. Possible options are 'logFC' or 'pvalue'. In case pvalue is selected, -log10(pval) will be multiplied by the sign of log2FC. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] dge_toptable gmt_file",
                     option_list = option_list, prog = "gsea",
                     description = "Perform log2FC-Preranked GSEA using the clusterProfiler R package. 
                     'dge_toptable' is the path of the DGE toptable from which genes are sorted according to log2FC. The DGE toptable must have the following structure:
                     .META: dge_toptable
                        1. gene names
                        2. baseMean (average expression)
                        3. log2FoldChange
                        4. pvalue
                        5. padj (FDR)
                        6... gene information (optional)
                     'gmt_file' is the path of a .gmt file containing pathways to evaluate. Each pathway must be reported in a different row, with genes separated by tabs.")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (dge_toptable and gmt_file)", call.=FALSE)
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


dge_toptable = arguments$args[1]
gmt_file = arguments$args[2]
rank = opt$ranking


# Importing dge toptable
toptable = read.delim(dge_toptable, h=T, row.names=1, check.names=F)
# Defining output name TODO
outname = gsub(".txt", "", gsub("deseq2_toptable", "gsea", dge_toptable))
collection = gsub(".gmt", "", gmt_file)
table_outname = paste(outname, collection, ".xlsx", sep="")
plots_outname = paste(outname, collection, ".pdf", sep="")


# Importing gmt
my_gmt = read.gmt(gmt_file)

# Ranked genelists (sorted according to logFC)
toptable = na.omit(toptable)
if(rank == "logFC"){
  my_rnk = toptable$log2FoldChange
  names(my_rnk) = rownames(toptable)
  my_rnk = sort(my_rnk, decreasing=T)
}
if(rank == "pvalue"){
  my_rnk = -log10(toptable$pvalue) * sign(toptable$log2FoldChange)
  names(my_rnk) = rownames(toptable)
  my_rnk = sort(my_rnk, decreasing=T)
}
if(rank!="logFC" & rank!="pvalue"){
  stop("The ranking statistics must be either 'logFC' or 'pvalue'")
}

# Running GSEA
my_gsea = GSEA(my_rnk, TERM2GENE=my_gmt, pvalueCutoff=0.5)
curated_gsea = as.data.frame(my_gsea)

if(nrow(curated_gsea)!=0){
  # Save results
  #xlsx::write.xlsx(curated_gsea, table_outname, sheetName=collection, col.names=T, row.names=F, append=F)
  openxlsx::write.xlsx(x=curated_gsea, file=table_outname)
  # Make plots
  plot_list = list()
  dim = ifelse(nrow(curated_gsea)<50, nrow(curated_gsea), 50)
  for(j in 1:dim){
    p = gseaplot2(my_gsea, geneSetID=j, title=my_gsea$Description[j])
    plot_list[[j]] = p
  }
  
  # Create pdf where each page is a separate plot
  pdf(plots_outname)
  for(j in 1:dim){print(plot_list[[j]])}
  dev.off()
} else {
  no_enrichment = data.frame(a="no_enrichment")
  #xlsx::write.xlsx(no_enrichment, table_outname, sheetName=collection, col.names=F, row.names=F, append=F)
  openxlsx::write.xlsx(x=no_enrichment, file=table_outname)
}




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()
