#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("SummarizedExperiment"))

### Script to produce violin plots of specific control genes 

option_list <- list(
  make_option(c("-a", "--annotations"), action="store", type="character", default="", help="The name of the categorical variable used to group samples into separate violin plots. It must be the name of the corresponding metadata column. [default \"%default\"]"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit."),
  make_option(c("-j", "--jitter"), action="store_true", default=FALSE, help="Whether to represent the single observations as jittered points.  [default \"%default\"]")
)

parser <- OptionParser(usage = "%prog [options] rna_object genes_list
                     'rna_object' is the path to a .rds object containing a Summarized Experiment with normalized expression data and samples metadata.
                     'genes_list' is the path to a tab-delimited file containing the list of control genes to evaluate. Each gene must be reported in a different row, with no header. Gene encoding must be the same used for the raw counts matrix.",
                     option_list = option_list, prog = "violin_plot",
                     description = "Produce violin plots on the normalized expression values (cpm) of specific control genes."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (rna_object and genes_list)", call.=FALSE)
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
genes_list = arguments$args[2]
annot = opt$annotations
jit = opt$jitter


# Importing rna_object
rna_exp = get(load(rna_object))
data = rna_exp@assays@data$cpm

# Importing control genes
control_genes = as.character(read.delim(genes_list, h=F)$V1)
control_genes = control_genes[control_genes%in%rownames(data)]
if(length(control_genes)==0){
  stop("None of the provided genes were retrieved in the expression matrix. Check consistency between the expression data and the supplied list of genes.")
}

# Filtering data table
data_ctrl = as.data.frame(t(data[control_genes,]))

if(annot == ""){
  # One single violin for all the samples
  pdf("violin_plots_control_genes.pdf")
  if(jit){
    for(i in colnames(data_ctrl)){
      p = ggplot(data_ctrl, aes(x="", y=data_ctrl[,i])) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_minimal()
      p = p + labs(title=i, y="CPM") + theme(title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))
      print(p)
    }
  } else {
    for(i in colnames(data_ctrl)){
      p = ggplot(data_ctrl, aes(x="", y=data_ctrl[,i])) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal()
      p = p + labs(title=i, y="CPM") + theme(title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))
      print(p)
    }
  }
  dev.off()
  
} else {
  if(!annot%in%colnames(colData(rna_exp))){
    stop("The selected variable does not appear in metadata. Please ensure the variable name matches exactly as it appears in the sample metadata.")
  }
  data_ctrl$group = colData(rna_exp)[,annot]
  if(sum(rownames(data_ctrl) != rownames(colData(rna_exp))) > 0){
    stop("Sample names inconsistency! Check the summarized experiment object.")
  }
  # One violin plot for each category of the selected variable
  pdf("violin_plots_control_genes.pdf")
  if(jit){
    for(i in colnames(data_ctrl)){
      p = ggplot(data_ctrl, aes(x=group, y=data_ctrl[,i], fill=group)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_minimal()
      p = p + labs(title=i, y="CPM") + theme(title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))
      print(p)
    }
  } else{
    for(i in colnames(data_ctrl)){
      p = ggplot(data_ctrl, aes(x=group, y=data_ctrl[,i], fill=group)) + geom_violin() + geom_boxplot(width=0.1) + theme_minimal()
      p = p + labs(title=i, y="CPM") + theme(title=element_text(size=20), axis.text=element_text(size=15), legend.text=element_text(size=15))
      print(p)
    }
  }
  dev.off()
}











########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()







