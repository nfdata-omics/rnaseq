#!/usr/bin/env Rscript

suppressMessages(library("optparse"))
suppressMessages(library("enrichR"))
suppressMessages(library("openxlsx"))

### Script to perform first line functional enrichment using enrichR 
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.05, help="The FDR cutoff to define significant genes. [default \"%default\"]"),
  make_option(c("-l", "--logfc"), action="store", type="double", default=0, help="The cutoff on log2FC absolute value to define significant genes (in combination with FDR). [default \"%default\"]")
)

parser<-OptionParser(usage = "%prog [options] dge_toptable",
                     option_list = option_list, prog = "ora_enrichr",
                     description = "Perform over-representation analysis using the enrichR R package. 
                     'dge_toptable' is the path of the DGE toptable from which significant genes are extracted. The DGE toptable must have the following structure:
                     .META: dge_toptable
                        1. gene names
                        2. baseMean
                        3. log2FoldChange
                        4. pvalue
                        5. padj (FDR)")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args)!=1) {
  stop("One argument must be supplied (dge_toptable)", call.=FALSE)
}

dge_toptable = arguments$args[1]
fdr = opt$FDR
logfc = opt$logfc


# Full list of pathway source
databases <- listEnrichrDbs()
# Selected databases to make the enrichment of (chosen from databases)
enrich.databases <- c("GO_Biological_Process_2023",
                      "GO_Cellular_Component_2023",
                      "GO_Molecular_Function_2023",
                      "Reactome_2022",
                      "KEGG_2021_Human",
                      "WikiPathways_2024_Human",
                      "BioCarta_2016",
                      "MSigDB_Hallmark_2020")
# Looking for more recent versions:
#databases$libraryName[grepl("Hallmark", databases$libraryName)]


# Importing dge toptable
toptable = read.delim(dge_toptable, h=T, row.names=1, check.names=F)
# Defining output name
outname = gsub("deseq2_toptable", "enrichr", dge_toptable)
outname = gsub(".txt", "", outname)

# Extracting significant genes using specified FDR and log2FC cutoffs
my_genes = rownames(toptable)[toptable$padj<fdr & abs(toptable$log2FoldChange)>logfc]
my_genes_up = rownames(toptable)[toptable$padj<fdr & toptable$log2FoldChange>logfc]
my_genes_down = rownames(toptable)[toptable$padj<fdr & toptable$log2FoldChange<(-logfc)]

# Performing enrichment and saving results
if(length(my_genes)>=5){
  my_enrichr = enrichR::enrichr(genes=my_genes, databases=enrich.databases)
  openxlsx::write.xlsx(x=my_enrichr, file=paste(outname,"_all.xlsx",sep=""))
}
if(length(my_genes_up)>=5){
  my_enrichr_up = enrichR::enrichr(genes=my_genes_up, databases=enrich.databases)
  openxlsx::write.xlsx(x=my_enrichr_up, file=paste(outname,"_up.xlsx",sep=""))
}
if(length(my_genes_down)>=5){
  my_enrichr_down = enrichR::enrichr(genes=my_genes_down, databases=enrich.databases)
  openxlsx::write.xlsx(x=my_enrichr_down, file=paste(outname,"_down.xlsx",sep=""))
}



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


