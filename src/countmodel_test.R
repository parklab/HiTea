#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Input files not specified. Usage: countmodel.R [Outprefix] [Workdir] [pval_cutoff] [TE_string] [RefMAPQ] [TEMAPQ]", call.=FALSE)
}
suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(MASS,quietly = T))
suppressPackageStartupMessages(library(data.table,quietly = T))
#args = c("test", "C:/d",0.05,"Alu,L1,SVA,HERVK",28,30)
#args = c("gm12878_20xsc", "C:/d/bgmodel",0.05,"Alu,L1,SVA,HERVK",28,30)

pvalcutoff = as.numeric(args[3])
tes = as.character(unlist(strsplit(args[4],",")))
refMAPQ = as.numeric(args[5])
TEMAPQ = as.numeric(args[6])
window=1000
RAMCUTOFF = 3
addPolyA=T ## Add PoyA read support


#####################################################################################################
a <- read.delim(paste0(args[2],"/",args[1],".temp.insertions.bed"),header=F,comment.char = "#")
names(a)[1:9] <- c("chr","start","end","name","score","strand","TE","id","cov")
a$rname = paste0(a$TE,"_",a$name)
a <- as(a,"GRanges")
a$cov <- a$cov*10

p <- read.delim(paste0(args[2],"/",args[1],"_RandomLocs.bed"),header=F,comment.char = "#")
names(p) <- c("chr","start","end","name","score","strand","cov")
p$rname <- paste0("random_",1:length(p$chr))
p <- as(p,"GRanges")
p <- resize(p,1)
p$cov <- p$cov*10
p <- p[p$cov>0]

########################## Genrate counts
## 1. Genreate separate count table for all potential insertions using all transposon supporting reads
## 2. Run background enrichment model for each class of transposon separately
## 3. Assign P-values to 1 if the count of supportng reads is <3 within 500bp of insertion site (typical assumption by WGS callers MELT and TEA)
getCounts <- function(gr,a1,refMAPQ=0,clipCov=T,window=1000,tot=F,ram=T,dircheck=T){ ## GenomicRanges and RangeObject 
  
  if(!"clip"%in%names(values((gr)))){
    stop("Clip location not provided in the GRanges object\n")
  }
  if(!"read_start"%in%names(values((gr)))){
    stop("read_start location not provided in the GRanges object\n")
  }
  if(!"read_end"%in%names(values((gr)))){
    stop("read_end location not provided in the GRanges object\n")
  }
  if(clipCov==F & ram==F){
    stop("Input parameters not specified correctly\n")
  }
  
  strand(a1) <- "+"
  detp <- gr[grep("DE,TP",gr$evi)]
  if(length(detp)>1){
    gr <- gr[-grep("DE,TP",gr$evi)]
  }
  gr <- gr[gr$refMAPQ>=refMAPQ]
  start(gr) <- gr$read_start
  end(gr) <- gr$read_end
  a1$totRam <- 0
  if(ram==T){
    if("TE"%in%names(values(gr))){
      gr <- gr[gr$TE!="PolyA"]
    }
    #if(dircheck==T){
    #  gr <- gr[strand(gr)!=as.character(gr$TEStrand)]
    #}
    a1$totRam <- countOverlaps( resize(a1,width(a1)+window,fix="center"),gr,ignore.strand=T )
  }
  if(length(detp)>1 & clipCov==T){
    start(detp) <- 1
    end(detp) <- 1
    start(detp) <- end(detp) <- detp$clip
    detp <- as(detp,"GRanges")
    if(tot==F){
      a1$totClip <- countOverlaps(a1,detp,ignore.strand=T)
    }else{
      a1$totClip <- countOverlaps(resize(a1,width(a1)+window,fix="center"),detp,ignore.strand=T)
    }
    a1$totRam <- a1$totRam+a1$totClip
  }
  
  return(a1$totRam)
}
model <- function(x){
  
  my_params <- "0,0"
  if(length(unique(x)) >2 & length(x) >=5){
    fit = tryCatch({
      eps <- sqrt(.Machine$double.eps)
      fitdistr(x,densfun = "negative binomial",lower=c(eps,eps))
    },
    error = function(e) {
      x <- rep(0,10)
      fitdistr(x,densfun = "negative binomial",lower=c(0.01,0.01))
    })
    my_params = paste0(fit$estimate[1],",",fit$estimate[2])
  }
  return(my_params)
}

cnt <- data.frame(id=p$rname)
cnta <- data.frame(id=a$rname)
GR <- GRanges()
for (f in c(tes,"PolyA","unmap")){
  cat(" counting non-conforming unsplit RPs in a window for ",f,"\n")
  f = gsub("/\\S*$","",f)
  load(paste0(args[2],"/",args[1],"_",f,".RData"))
  #if(f != "unmap"){
    cnt <- cbind(cnt,getCounts(gr,p,refMAPQ=refMAPQ,clipCov=F,window=4000,tot=F,dircheck = F))
    cnta <- cbind(cnta,getCounts(gr,a,refMAPQ=refMAPQ,clipCov=F,window=4000,tot=F,dircheck = F))
  #}
  start(gr) <- 1  ## to omit the error message
  end(gr) <- 1  ## to omit the error message
  start(gr) <- end(gr) <- gr$clip
  GR <- c(GR,gr[grep("DE,TP",gr$evi),"id"])
  rm(gr,f)
}
cat("\n")
names(cnta)[2:ncol(cnta)] <- names(cnt)[2:ncol(cnt)] <- c(tes,"PolyA","unmap")
a$disc <- 0
p$disc <- 0
#for(f in tes){
#  a[grep(f,a$rname)]$disc <- cnta[grep(f,cnta$id),f]
#  p[grep(f,p$rname)]$disc <- cnt[grep(f,cnt$id),f]
#}
p$disc <- rowSums(cnt[,2:ncol(cnt)])
a$disc <- rowSums(cnta[,2:ncol(cnta)])
GR <- data.table(as.data.frame(GR))
setkeyv(GR, c('id'))
GR <- subset(unique(GR))
GR <- as.data.frame(GR)
GR <- as(GR,"GRanges")
a$fzcnt <- countOverlaps(resize(a,width(a)+100,"center"),GR,ignore.strand=T )
rm(cnt,cnta,GR)

print.figure=T
if(print.figure==T){
  pdf(file= paste0(args[2],"/",args[1],"_bgModeling.pdf"),width = 12,height = 16)
  par(mfrow=c(length(tes),3))
}
res <- data.frame(id=a$rname)
plotdf <- list()
bw=1

for (f in tes){
  f = gsub("/\\S*$","",f)
  load(paste0(args[2],"/",args[1],"_",f,".RData"))
  gr$TE <- paste0(f)
  
  if(addPolyA==T & file.exists(paste0(args[2],"/",args[1],"_","PolyA",".RData")) ){  
    gr1 <- gr
    rm(gr)
    load(paste0(args[2],"/",args[1],"_","PolyA",".RData"))
    gr$TE <- "PolyA"
    gr <- c(gr,gr1)
    rm(gr1)
    cat(" creating background model and estimating enrichment for ", f," (w PolyA)\t")
  }else{
    cat(" creating background model and estimating enrichment for ", f," (wo PolyA)\t")
  }
  
  ctr <- data.frame(id=p$rname,x=getCounts(gr,p,refMAPQ=refMAPQ,clipCov=T,window=window,tot=F,ram=T),y=round(p$disc))
  ex <- data.frame(id=a$rname, x=getCounts(gr,a,refMAPQ=refMAPQ,clipCov=T,window=window,tot=F,ram=T),y=round(a$disc))
  
  #if(nrow(ctr[ctr$y<ctr$x,])>0 | nrow(ex[ex$y<ex$x,])>0 ){
  #  stop("error encountered while counting the non-conforming reads\n")
  #}
  ctr$ram <- getCounts(gr,p,refMAPQ=refMAPQ,clipCov=F,window=window,tot=F,ram=T)
  ex$ram <- getCounts(gr,a,refMAPQ=refMAPQ,clipCov=F,window=window,tot=F,ram=T)
  
  ctr$y <- floor(ctr$y/bw)
  ex$y <- floor(ex$y/bw)
  ctr <- data.table(ctr)
  ctr <- ctr[,parameters:= model(x),by=list(y)]
  ctr <- ctr[,sd:= sd(x),by=list(y)]
  ctr <- ctr[,mean:= mean(x),by=list(y)]
  ctr <- as.data.frame(ctr)
  bg <- unique(ctr[,c("y","parameters","sd","mean")])
  
  z <- unique(ex$y)
  z <- z[!z%in%bg$y]
  if(length(z)>0){
    bg <- rbind(bg, data.frame(y=z,parameters="0,0",sd=NA,mean=NA))
  }
  bg <- bg[order(bg$mean,decreasing = F),]
  
  ## plot expected Poisson
  pl <- data.frame(x=(ctr$y+ctr$x),y=ctr$x,dataframe="ctr")
  pl <- rbind(pl, data.frame(x=bg$mean,y=bg$sd,dataframe="bg"))
  
  if(print.figure==T){
    smoothScatter((ctr$y+ctr$x),(ctr$x),main=paste0(f),xlab="Total coverage",ylab="Hi-C disc. read coverage",cex.lab=1.5,cex.axis=1.5,cex.main=3)
    loess.fit = tryCatch( {loess.smooth((ctr$y+ctr$x),(ctr$x)) },
                          error = function(e) { 
                            l <- list()
                            return(l)}
    )
    if(length(loess.fit)>0){
      lines(loess.fit,col="red")
    }
    plot(bg$mean,bg$sd,main=paste0(f),xlab="mean",ylab="sd",cex=1.2,cex.lab=1.5,cex.axis=1.5,cex.main=3,pch=16,col=rgb(1,0,0,0.2))
    lines(x = bg$mean,y=sqrt(bg$mean),col="black")
    legend("bottomright",c("poisson expectation"),col=c("black"),lwd=2,border = NULL)
  }
  
  ## Extrapolate the parameters 
  bg <- bg[order(bg$y,decreasing = F),]
  ol <- "0,0"
  for(i in 1:length(bg$y)){
    if(as.character(bg$parameters[i]) != "0,0"){
      ol <- as.character(bg$parameters[i])
    }
    if(as.character(bg$parameters[i]) == "0,0" ){
      bg$parameters[i] <- ol
    }
  }
  
  ## merge bg model parameters with ex
  ex <- merge(ex,bg[,1:2],by="y",all.x=T)
  ex$p1 <- as.numeric(gsub(",\\S*","",ex$parameters))
  ex$p2 <- as.numeric(gsub("\\S*?,","",ex$parameters))
  ex$p.val <- apply(ex[,c("x","p1","p2")],1,function(x){
    if( (x[2]==0 & x[3]==0) ){
      return(1)
    }else{
      return(pnbinom(x[1],size =x[2],mu = x[3],lower.tail = F ))
    }
  })
  
  ctrx <- merge(ctr[,c("x","y")],bg[,c("y","parameters")],by="y",all.x=T)
  ctrx$p1 <- as.numeric(gsub(",\\S*","",ctrx$parameters))
  ctrx$p2 <- as.numeric(gsub("\\S*?,","",ctrx$parameters))
  ctrx$p.val <- apply(ctrx[,c("x","p1","p2")],1,function(x){
    if( (x[2]==0 & x[3]==0) ){
      return(1)
    }else{
      return(pnbinom(x[1],size =x[2],mu = x[3],lower.tail = F ))
    }
  })
  if(print.figure==T){
    hist(ctrx$p.val,breaks=25,xlab="P.value",ylab = "Frequency",main=paste0(f," (p-values)"),cex.lab=1.5,cex.axis=1.5,cex.main=3,xlim=c(0,1))
    #hist(ex$p.val,breaks=25,xlab="P.value",ylab = "Frequency",main=paste0(f," (p-values)"),cex.lab=1.5,cex.axis=1.5,cex.main=3,xlim=c(0,1))
  }
  pl <-rbind(pl,data.frame(x=ex$p.val,y=1,dataframe="ex"))
  plotdf[[paste0(f)]] <- pl  
  
  q1 <- ex[ex$p.val<pvalcutoff,]
  q1 <- q1[grep(f,q1$id),]
  cat(": ",length(q1$id),"putative cases (pval<",pvalcutoff,")\n")
  ex$p.val <- ifelse(ex$x<RAMCUTOFF,1,ex$p.val)
  
  ex <- ex[match(res$id,ex$id),]
  ex <- ex[,c("p.val","ram")]
  names(ex) <- c("p.val","reads")
  res <- cbind(res,ex)
  rm(gr,f,ctr,ex,bg,z,pl,loess.fit,q1,i,ctrx,ol)
}

if(print.figure==T){
  garbage <- dev.off()
}

res$propable.TE <- apply(res,1,function(y){
  #y <- unname(as.vector(y))
  id = as.character(unlist(y[1]))
  id = gsub("_\\d*","",id)
  x <- as.numeric(y[seq(2,length(y)-1,2)]) ## pvalues
  z <- as.numeric(y[seq(3,length(y),2)]) ## RAM counts
  names(z) <- names(x) <- 1:length(x)
  x <- x[x<1]
  x <- sort(x,decreasing = F)
  z <- z[names(z)%in%names(x)] 
  z <- z[match(names(x),names(z))]
  x <- formatC(x, format = "e", digits = 2)
  
  w <- "BG="
  wa <- "Rest_BG="
  for(i in names(x)){
    if(id == as.character(tes[as.numeric(i)])){
      w <- paste0(w,paste0(tes[as.numeric(i)],",",x[i],",",z[i],";"))
    }else{
      wa <- paste0(wa,paste0(tes[as.numeric(i)],",",x[i],",",z[i],";"))
    }
  }
  
  if((length(x)>1 & id =="*") |length(x)==0){
    w <- "BG=;"
    wa <- "Rest_BG=;"
  }
  if((length(x)==1) & id == "*"){
    w <- wa
    wa <- "Rest_BG=;"
    w <- gsub("Rest_","",w)
  }
  if(wa =="Rest_BG="){
    wa <- "Rest_BG=;"
  }
  if(w == "BG="){
    w <- "BG=;"
  }
  o <- paste0(w,wa)
  o
}) 
a$bg.model = res$propable.TE
a$disc <- NULL

# test the performance
if(T){
  myf<- function(pacb){
    h <- distanceToNearest(pacb) %>% as.data.frame()
    h$ref <- pacb[h$queryHits]$TE
    h$test <- pacb[h$subjectHits]$TE
    h <- h[h$distance<100,]
    h <- h[h$ref==h$test,]
    h$queryHits %>% unique
    h <- h[duplicated(t(apply(h[,1:2],1,sort))),]
    if(nrow(h)>0){
      pacb <- pacb[-h$queryHits]
    }
    pacb
  }
  pacb = rbind(read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_pacbio_LINE1.txt",header=F),
               read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_pacbio_Alu.txt",header=F),
               read.delim("A:/work/RT in HiC/xTEA_simon/xTEA_pacbio_SVA.txt",header=F))
  pacb <- with(pacb,GRanges(V1,IRanges(V2,V2),"*",TE=V3))
  pacb$TE <- gsub("LINE1","L1",pacb$TE)
  alu <- myf(unique(pacb[pacb$TE=="Alu"]))
  sva <- myf(unique(pacb[pacb$TE=="SVA"]))
  l1 <- myf(unique(pacb[pacb$TE=="L1"]))
  pacb<- c(alu,l1,sva)
  table(pacb$TE)
  rm(alu,sva,l1)
  
  a$pval <- gsub(",.*$","",gsub("^.*?,","",a$bg.model) ) %>% as.numeric
  sum(is.na(a$pval))
  a$pval <- ifelse(is.na(a$pval),1,a$pval)
  a$alu <- countOverlaps(resize(a,101,"center"),pacb[pacb$TE=="Alu"],ignore.strand=T)
  a$l1 <- countOverlaps(resize(a,101,"center"),pacb[pacb$TE=="L1"],ignore.strand=T)
  a$sva <- countOverlaps(resize(a,101,"center"),pacb[pacb$TE=="SVA"],ignore.strand=T)
  
  #a$pval <- p.adjust(a$pval,method = "fdr")
  
  cat(sum(a$TE=="Alu" & a$pval<0.05 & a$alu>0),"\n")
  cat(sum(a$TE=="L1" & a$pval<0.05 & a$l1>0),"\n")
  cat(sum(a$TE=="SVA" & a$pval<0.05 & a$sva>0),"\n\n")
  
  cat(sum(a$TE=="Alu" & a$pval<0.05 & a$alu==0),"\n")
  cat(sum(a$TE=="L1" & a$pval<0.05 & a$l1==0),"\n")
  cat(sum(a$TE=="SVA" & a$pval<0.05 & a$sva==0),"\n\n")
  
}


a = as.data.frame(a)
save(plotdf,file= paste0(args[2],"/",args[1],"_bgModeling.RData"))
write.table(a,file=paste0(args[2],"/",args[1],".temp.bgModeledInsertions.txt"),quote = F,sep = "\t",row.names = F,col.names = F)


