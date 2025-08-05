#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("enrichplot"))
suppressMessages(library("msigdbr"))
suppressMessages(library("openxlsx"))

### Script to perform first line functional enrichment using over-representation analysis from clusterProfiler
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.05, help="The FDR cutoff to define significant genes. [default \"%default\"]"),
  make_option(c("-l", "--logfc"), action="store", type="double", default=0, help="The cutoff on log2FC absolute value to define significant genes (in combination with FDR). [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] dge_toptable gmt_file",
                     option_list = option_list, prog = "ora_clusterProfiler",
                     description = "
                     Perform over-representation analysis using the clusterProfiler R package. 
                     'dge_toptable' is the path of the DGE toptable from which significant genes are extracted. The DGE toptable must have the following structure:
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
fdr = opt$FDR
logfc = opt$logfc


# Importing dge toptable
toptable = read.delim(dge_toptable, skip=1, h=T, row.names=1, check.names=F)
# Defining output name
outname = gsub(".txt", "", gsub("deseq2_toptable", "ora_CP", dge_toptable))
collection = gsub(".gmt", "", gmt_file)
#table_outname = paste(outname, collection, "xlsx", sep=".")
#plots_outname = paste(outname, collection, "pdf", sep=".")

# Extracting significant genes using specified FDR and log2FC cutoffs
background_genes = rownames(toptable)
toptable = toptable[!is.na(toptable$padj),]
my_genes = rownames(toptable)[toptable$padj<fdr & abs(toptable$log2FoldChange)>logfc]
my_genes_up = rownames(toptable)[toptable$padj<fdr & toptable$log2FoldChange>logfc]
my_genes_down = rownames(toptable)[toptable$padj<fdr & toptable$log2FoldChange<(-logfc)]

# Importing gmt
my_gmt = read.gmt(gmt_file)

# Enrichment
if(length(my_genes)>=5){
  if(length(my_genes)>10000){
    my_genes = my_genes[1:10000]
    print("The list of differentially expressed genes was longer than 10.000 genes, therefore it was truncated at 10.000 genes. 
            Consider to apply more stringent cutoffs on FDR and/or log2FC")
  }
  my_enrich = tryCatch({ 
    enricher(gene=my_genes, universe=background_genes, maxGSSize=1000, TERM2GENE=my_gmt)},
    # BgRatio = M/N, with M = #genes in the selected pathway that are also in the background; with N = #unique genes in all the gmt (all pathways) that are also in the background
    # GeneRatio = M/N, with M = #my input genes that are also in the selected pathway; with N = #my input genes that are also in all the gmt (all pathways)
    error = function(e) {message(e$message)
  })
  if (!is.null(my_enrich)) {
    write.table(x=my_enrich@result, file=paste(outname, ".", collection, ".all.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  } else {
    no_enrichment = data.frame(empty="No enrichment was found.")
    write.table(x=no_enrichment, file=paste(outname, ".", collection, ".all.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
} else {
  no_enrichment = data.frame(empty=paste("Enrichment analysis was not performed. Less than 5 significant genes were retrieved using the cutoff FDR<", fdr, ", abs_log2FC>", logfc, ".", sep=""))
  write.table(x=no_enrichment, file=paste(outname, ".", collection, ".all.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

if(length(my_genes_up)>=5){
  if(length(my_genes_up)>10000){
    my_genes_up = my_genes_up[1:10000]
    print("The list of up-regulated genes was longer than 10.000 genes, therefore it was truncated at 10.000 genes. 
            Consider to apply more stringent cutoffs on FDR and/or log2FC")
  }
  my_enrich_up = tryCatch({ 
    enricher(gene=my_genes_up, universe=background_genes, maxGSSize=1000, TERM2GENE=my_gmt)},
    error = function(e) {message(e$message)
  })
  if(!is.null(my_enrich_up)){
    write.table(x=my_enrich_up@result, file=paste(outname, ".", collection, ".up.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  } else {
    no_enrichment = data.frame(empty="No enrichment was found.")
    write.table(x=no_enrichment, file=paste(outname, ".", collection, ".up.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
} else {
  no_enrichment = data.frame(empty=paste("Enrichment analysis was not performed. Less than 5 significant genes were retrieved using the cutoff FDR<", fdr, ", abs_log2FC>", logfc, ".", sep=""))
  write.table(x=no_enrichment, file=paste(outname, ".", collection, ".up.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}

if(length(my_genes_down)>=5){
  if(length(my_genes_down)>10000){
    my_genes_down = my_genes_down[1:10000]
    print("The list of down-regulated genes was longer than 10.000 genes, therefore it was truncated at 10.000 genes. 
            Consider to apply more stringent cutoffs on FDR and/or log2FC")
  }
  my_enrich_down = tryCatch({
    enricher(gene=my_genes_down, universe=background_genes, maxGSSize=1000, TERM2GENE=my_gmt)},
    error = function(e) {message(e$message)
  })  
  if(!is.null(my_enrich_down)){
    write.table(x=my_enrich_down@result, file=paste(outname, ".", collection, ".down.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  } else {
    no_enrichment = data.frame(empty="No enrichment was found.")
    write.table(x=no_enrichment, file=paste(outname, ".", collection, ".down.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
  }
} else {
  no_enrichment = data.frame(empty=paste("Enrichment analysis was not performed. Less than 5 significant genes were retrieved using the cutoff FDR<", fdr, ", abs_log2FC>", logfc, ".", sep=""))
  write.table(x=no_enrichment, file=paste(outname, ".", collection, ".down.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F)
}





########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


