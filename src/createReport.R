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
rmarkdown::render('HiTEA_Report.Rmd',params = list(dir=args[2],outprefix=args[1]),output_dir = args[2])
