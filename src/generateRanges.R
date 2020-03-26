#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Input supporting read file not provided", call.=FALSE)
} else if(length(args)==1){
  args[2] = paste0(args[1])
}
args[2] = paste0(args[2],".RData")

suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(data.table,quietly = T))
##(1)chr,(2)start,(3)end,(4)id,(5)strand,(6)evi,(7)clip, (8)refMapqQ, (9)TEMapScore, (10) TE_strand (11) TE_mapPos (12) TE
gr <- read.delim(args[1],header=F)
gr <- data.table(gr)
setkeyv(gr, c('V4','V6'))
gr <- subset(unique(gr))
gr <- as.data.frame(gr)
gr <- with(gr,GRanges(V1,IRanges(V2,V3),V5,evi=V6,clip=V7,read_start=V2,read_end=V3,refMAPQ=V8,TEMAPQ=V9,TEStrand=V10,TEMap=V11,id=V4))
gr <- resize(gr,1,"start")
gr <- sortSeqlevels(gr)
gr <- sort(gr)
save(gr,file=args[2])
