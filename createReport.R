#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Input files not specified. Usage: createReport.R [Outprefix] [Workdir]", call.=FALSE)
}
lib_path="/n/data1/hms/dbmi/park/Dhawal/R"
suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(data.table,quietly = T))
suppressPackageStartupMessages(library(knitr,quietly = T,lib.loc = lib_path))
suppressPackageStartupMessages(library(rmarkdown,quietly = T,lib.loc = lib_path))
suppressPackageStartupMessages(library(ggplot2,quietly = T,lib.loc = lib_path))
suppressPackageStartupMessages(library(EnrichedHeatmap,quietly = T,lib.loc = lib_path))
suppressPackageStartupMessages(library(circlize,quietly = T,lib.loc = lib_path))
#suppressPackageStartupMessages(library(kableExtra,quietly = T,lib.loc = lib_path))
#suppressPackageStartupMessages(library(DT,quietly = T,lib.loc = lib_path))

## wrapper function

rmarkdown::render('markdownfile.rmd',params = list(dir=args[2],outprefix=args[1]),output_dir = args[2])

#install.packages('kableExtra',lib = "/n/data1/hms/dbmi/park/Dhawal/R3.3",repos = "http://cran.us.r-project.org",dependencies = T)
#install.packages(c('DT',"ggplot2","data.table","knitr","rmarkdown","dplyr","circlize","reshape2","MASS"),lib = "/n/data1/hms/dbmi/park/Dhawal/R3.3",repos = "http://cran.us.r-project.org",dependencies = T)
#rmarkdown::render('markdownfile.rmd',params = list(dir="C:/d/report",outprefix="gm12878_20x"),output_dir = "C:/d/report")
