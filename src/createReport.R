#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Input files not specified. Usage: createReport.R [Outprefix] [Workdir]", call.=FALSE)
}
suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(data.table,quietly = T))
suppressPackageStartupMessages(library(knitr,quietly = T))
suppressPackageStartupMessages(library(rmarkdown,quietly = T))
suppressPackageStartupMessages(library(ggplot2,quietly = T))
suppressPackageStartupMessages(library(EnrichedHeatmap,quietly = T))
suppressPackageStartupMessages(library(circlize,quietly = T))
#suppressPackageStartupMessages(library(kableExtra,quietly = T,lib.loc = lib_path))
#suppressPackageStartupMessages(library(DT,quietly = T,lib.loc = lib_path))

rmarkdown::render('HiTEA_Report.Rmd',params = list(dir=args[2],outprefix=args[1]),output_dir = args[2])
#rmarkdown::render('A:/work/scripts/TE_insertions/final/v1.3/markdownfile.Rmd',params = list(dir="C:/d/temp/plots",outprefix="gm12878_20xsb"),output_dir = "C:/d/temp/plots")




