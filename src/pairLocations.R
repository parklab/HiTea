#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Input cluster file not specified", call.=FALSE)
} else if(length(args)==1){
  args[2] = "temp.cluster.template.txt"
}

suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(data.table,quietly = T))

#d <- read.delim("A:/work/temp/gmH1.temp.clusters.txt",header=T)
#d <- with(d,GRanges(V1,IRanges(V4,V5),V3,V2))
d = read.delim(args[1],header=T,comment.char = "#")   ## chr, start,end,strand,cliploc,side
d = with(d,GRanges(chr,IRanges(start,end),strand,clipcoord,side,reads,score))
r = reduce(d,ignore.strand=T)
values(r) = cbind(values(r),data.frame(id=1:length(r)))
olap = findOverlaps(d,r,ignore.strand=T,select="all")
res = cbind(id=r[subjectHits(olap)]$id,as.data.frame(d[queryHits(olap)]))

## Approach:
### 1. get cluster with maximum reads
### 2. select 1st and 2nd location based on maximum score
### 3. For a side reads, use minimum cliplocation
### 4. for b side reads, use maximum cliplocation
res = data.table(res)
res = res[,a:=max(reads),by=list(id,seqnames,side)]
res = res[res$reads==res$a,]
res = res[,b:=max(score),by=list(id,seqnames,side,a)]
res = res[res$score==res$b,]
resa = res[res$side=="a",]
resb = res[res$side=="b",]
resa = resa[,coord:=min(clipcoord),by=list(id,seqnames,a,b)]
resb = resb[,coord:=max(clipcoord),by=list(id,seqnames,a,b)]

res = rbind(resa,resb)
res = as.data.frame(res)
res = unique(res[,c("id","seqnames","coord","side","a","b")])
write.table(res,file=args[2],sep="\t",quote = F,row.names = F,col.names = F)

###############################################################################################