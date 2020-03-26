#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Input files not specified. Usage: makeMat.R [Outprefix] [Workdir] [TE string]", call.=FALSE)
}
suppressPackageStartupMessages(library(GenomicRanges,quietly = T))

m.spl <- c()
m.unspl <- c()
WINDOW=10010
bw=10
tes <- as.character(unlist(strsplit(args[3],",")))

a <- read.delim(paste0(args[2],"/",args[1],".candidate.insertions.bed"),header=F,comment.char = "#")
names(a) <- c("chr","start","end","id","score","strand","TE","Evi","Descr","remark")
rownames(a) <- a$id
a <- as(a,"GRanges")

for (f in tes){
  load(paste0(args[2],"/",args[1],"_",f,".RData"))
  cat("generating matrix for ", f,"\n")
  gr$TE <- paste0(f)
  if(T){  ## Add PoyA read support
    gr1 <- gr
    rm(gr)
    if(file.exists(paste0(args[2],"/",args[1],"_","PolyA",".RData"))){
      load(paste0(args[2],"/",args[1],"_","PolyA",".RData"))
      gr$TE <- "PolyA"
      gr1 <- c(gr,gr1)
      rm(gr)
    }
    gr <- gr1
    rm(gr1)
  }
  
  
  gr <- resize(gr,1)
  gr.us <- gr[grep("IE,FP|IE,TP|DE,FP",gr$evi)]
  gr.us <- gr.us[gr.us$TE==f]
  #gr.us <- gr.us[gr.us$refMAPQ>=refMAPQ & gr.us$TEMAPQ>=TEMAPQ]
  gr.spl <- gr[grep("DE,TP",gr$evi)]
  
  a1 <- a[a$TE==f]
  ind <- unlist(tile(resize(a1,WINDOW,"center"),width = bw ))
  vec <- countOverlaps(ind,gr.spl,ignore.strand=T)
  x <- as.data.frame(matrix(vec,nrow = length(a1),ncol = round(WINDOW/bw), byrow = T) )
  rownames(x) <- a1$id
  if(nrow(x)>0){
    m.spl <- rbind(m.spl, x)
  }
  
  vec <- countOverlaps(ind,gr.us,ignore.strand=T)
  y <- as.data.frame(matrix(vec,nrow = length(a1),ncol = round(WINDOW/bw),byrow = T) )
  rownames(y) <- a1$id
  if(nrow(y)>0){
    m.unspl <- rbind(m.unspl, y)
  }
  rm(gr,gr.us,gr.spl,x,y,vec,a1)
}
m.spl <- m.spl[match(a$id,rownames(m.spl)),]
m.unspl <- m.unspl[match(a$id,rownames(m.unspl)),]
rownames(m.spl) <- rownames(m.unspl) <- a$id

save(m.spl,m.unspl,file=paste0(args[2],"/",args[1],".CovPlotMatrix.RData"))
