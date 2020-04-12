#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Input files not specified. Usage: createReport.R [Outprefix] [Workdir] [BASEDIR]", call.=FALSE)
}
suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(data.table,quietly = T))
suppressPackageStartupMessages(library(knitr,quietly = T))
suppressPackageStartupMessages(library(rmarkdown,quietly = T))
suppressPackageStartupMessages(library(ggplot2,quietly = T))
suppressPackageStartupMessages(library(EnrichedHeatmap,quietly = T))
suppressPackageStartupMessages(library(circlize,quietly = T))
suppressPackageStartupMessages(library(kableExtra,quietly = T))
suppressPackageStartupMessages(library(DT,quietly = T))

if(args[3] ==""){
 rmarkdown::render(paste0('src/HiTEA_Report.Rmd'),params = list(dir=args[2],outprefix=args[1]),output_dir = args[2])
}else{
 rmarkdown::render(paste0(args[3],'/src/HiTEA_Report.Rmd'),params = list(dir=args[2],outprefix=args[1]),output_dir = args[2])
}

#rmarkdown::render('src/HiTEA_Report.Rmd',params = list(dir='/usr/local/bin/test_out',outprefix='test'),output_dir = '/usr/local/bin/tset_out')
#rmarkdown::render('A:/work/scripts/TE_insertions/final/v1.3/markdownfile.Rmd',params = list(dir="C:/d/temp/plots",outprefix="gm12878_20xsb"),output_dir = "C:/d/temp/plots")




