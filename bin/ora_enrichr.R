#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("enrichR"))
suppressMessages(library("openxlsx"))

### Script to perform first line functional enrichment using enrichR 
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.05, help="The FDR cutoff to define significant genes. [default \"%default\"]"),
  make_option(c("-l", "--logfc"), action="store", type="double", default=0, help="The cutoff on log2FC absolute value to define significant genes (in combination with FDR). [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] dge_toptable",
                     option_list = option_list, prog = "ora_enrichr",
                     description = "
                     Perform over-representation analysis using the enrichR R package. 
                     'dge_toptable' is the path of the DGE toptable from which significant genes are extracted. The DGE toptable must have the following structure:
                     .META: dge_toptable
                        1. gene names
                        2. baseMean (average expression)
                        3. log2FoldChange
                        4. pvalue
                        5. padj (FDR)")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=1 & !version) {
  stop("One argument must be supplied (dge_toptable)", call.=FALSE)
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
fdr = opt$FDR
logfc = opt$logfc


# Full list of pathway source
#databases <- listEnrichrDbs()
# Selected databases to make the enrichment of (chosen from databases)
enrich.databases <- c("GO_Biological_Process_2023",
                      "GO_Cellular_Component_2023",
                      "GO_Molecular_Function_2023",
                      "Reactome_2022",
                      "KEGG_2021_Human",
                  #    "WikiPathways_2024_Human",
                      "BioCarta_2016",
                      "MSigDB_Hallmark_2020")
# Looking for more recent versions:
#databases$libraryName[grepl("Hallmark", databases$libraryName)]


# Importing dge toptable
toptable = read.delim(dge_toptable, h=T, row.names=1, check.names=F)
# Defining output name
outname = gsub(".txt", "", gsub("deseq2_toptable", "enrichr", dge_toptable))
summary_name = gsub("deseq2_toptable", "enrichr_summary", dge_toptable)

# Extracting significant genes using specified FDR and log2FC cutoffs
toptable = toptable[!is.na(toptable$padj),]
my_genes = rownames(toptable)[toptable$padj<fdr & abs(toptable$log2FoldChange)>logfc]
my_genes_up = rownames(toptable)[toptable$padj<fdr & toptable$log2FoldChange>logfc]
my_genes_down = rownames(toptable)[toptable$padj<fdr & toptable$log2FoldChange<(-logfc)]

# Performing enrichment and saving results
signif_enrich_all = 0;
signif_enrich_up = 0;
signif_enrich_down = 0;
if(length(my_genes)>=5){
  if(length(my_genes)>10000){
    my_genes = my_genes[1:10000]
    print("The list of differentially expressed genes was longer than 10.000 genes, therefore it was truncated at 10.000 genes. 
            Consider to apply more stringent cutoffs on FDR and/or log2FC")
  }
  my_enrichr = enrichR::enrichr(genes=my_genes, databases=enrich.databases)
  openxlsx::write.xlsx(x=my_enrichr, file=paste(outname,"_all.xlsx",sep=""))
  signif_enrich_all = sum(unlist(lapply(my_enrichr, function(x) return(sum(x$Adjusted.P.value<0.05)))))
}
if(length(my_genes_up)>=5){
  if(length(my_genes_up)>10000){
    my_genes_up = my_genes_up[1:10000]
    print("The list of up-regulated genes was longer than 10.000 genes, therefore it was truncated at 10.000 genes. 
            Consider to apply more stringent cutoffs on FDR and/or log2FC")
  }
  my_enrichr_up = enrichR::enrichr(genes=my_genes_up, databases=enrich.databases)
  openxlsx::write.xlsx(x=my_enrichr_up, file=paste(outname,"_up.xlsx",sep=""))
  signif_enrich_up = sum(unlist(lapply(my_enrichr_up, function(x) return(sum(x$Adjusted.P.value<0.05)))))
}
if(length(my_genes_down)>=5){
  if(length(my_genes_down)>10000){
    my_genes_down = my_genes_down[1:10000]
    print("The list of up-regulated genes was longer than 10.000 genes, therefore it was truncated at 10.000 genes. 
            Consider to apply more stringent cutoffs on FDR and/or log2FC")
  }
  my_enrichr_down = enrichR::enrichr(genes=my_genes_down, databases=enrich.databases)
  openxlsx::write.xlsx(x=my_enrichr_down, file=paste(outname,"_down.xlsx",sep=""))
  signif_enrich_down = sum(unlist(lapply(my_enrichr_down, function(x) return(sum(x$Adjusted.P.value<0.05)))))
}

enrich_summary = paste("There are ", signif_enrich_up, " pathways significantly enriched with up-regulated genes.
                There are ", signif_enrich_down, " pathways significantly enriched with down-regulated genes.
                There are ", signif_enrich_all, " pathways significantly enriched with overall deregulated genes.", sep="")
writeLines(enrich_summary, con=summary_name)




########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()


