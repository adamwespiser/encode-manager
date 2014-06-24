library(hexbin)
library(ggplot2)
library(plyr)
library(reshape)
library(hexbin)
library(RColorBrewer)
library(reshape)
library(ggdendro)
x <- "yea this works..."
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)

withinRange3 <- function(x,lower,upper){
  if(x >= lower && x <= upper ){
    TRUE
  }
  else{
    FALSE
  }
  
}


entropy <- function(...){
  x = c(...)
  x = x / sum(x)
  y = log2(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}

entropyVec <- function(x,logFn=log2){
  x = x / sum(x)
  y = logFn(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}

JS <- function(x,y,logFn=log2){
  x = x / sum(x)
  y = y / sum(y)
  a = (x + y)/2
  entropyVec(a ,logFn)- ((entropyVec(x,logFn)+entropyVec(y,logFn))/2) -> z
#  entropyVec(a )- ((entropyVec(x)+entropyVec(y))/2) -> z
  z}

JSsp <- function(e1,e2,logFn=log2){
  1 - sqrt(JS(e1,e2,logFn))
}

calcTissSpec <- function(x,cols,logFn=log2){
  x <- sapply(x,as.numeric)
  profile <- diag(length(cols))
  specEn <- sapply(1:length(cols),function(y)JSsp(x,profile[y,],logFn))
  # cols[which(specEn == max(specEn))]
  max(specEn)
  #specEn
}

calcTissSpecVector <- function(x,logFn=log2){
  x <- sapply(x,as.numeric)
  profile <- diag(length(x))
  specEn <- sapply(seq_along(x),function(y)JSsp(x,profile[y,],logFn))
  # cols[which(specEn == max(specEn))]
  max(specEn)
  #specEn
}

#this needs to be generalized...
tissSpec <- function(...){
  x <- unlist(...)
  #  print(x)
  profile = list(c(1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1))
  specEn = c()
  for (i in 1:5){
    specEn[i] = JSsp(x,profile[[i]])
  }
  #print(specEn)
  c("HELAS3","GM12878","H1HESC","HEPG2","K562")[which(specEn== max(specEn))]
  max(specEn)
}


logNormal <- function(X,c=1){
  normX <- min(X) + 1
  log(normX)
}


maxExpr <- function(x,cols){
  ifelse( sum(x) > 0.00000000001,
          cols[which(x== max(x))],
          "none"
  )
}

threeTimesRestSpecifity <- function(expr.vec,cols){

  ifelse( sort(expr.vec,decreasing=TRUE)[1] > (3* sort(expr.vec,decreasing=TRUE)[2]),
          cols[which(expr.vec == max(expr.vec))],
          "none")}

nTimesRestSpecifity <- function(expr.vec,cols){
  first <-  sort(expr.vec,decreasing=TRUE)[1]
  second <-  sort(expr.vec,decreasing=TRUE)[2]
  n <- Inf # case that first >0, second==0
  if (first == 0 && second == 0){
     n <- 0
  }
  if (first > 0 && second> 0){
     n <- first/second
  }
  n
}



applyAndAppendFnDf <-function(df=df,cols=cols,FUN=fn,newCol=newCol){
  df[newCol] <- apply(df[as.vector(cols)],1,function(...){x<-c(...);FUN(x)})
  df
}

applyAndAppendFnDfSecondArg <-function(df=df,cols=cols,FUN=fn,newCol=newCol,arg2=arg2){
  df[newCol] <- apply(df[as.vector(cols)],1,function(...){x<-c(...);FUN(x,arg2)})
  df
}

#split vec 
exprLabels <- function(x){
 # y<-cut(x, quantile(x[which( x> 0)],(0:3)/3),labels=c("low","mid","high"))
  output.vec <- rep("none",length(x))
  expr.labels.wNA <-  as.character(cut(x, quantile(x[which( x> 0)],(0:3)/3),labels=c("low","mid","high")))
  expr.labels.index <- !is.na(expr.labels.wNA)
  output.vec[expr.labels.index] <- expr.labels.wNA[expr.labels.index]
  factor(output.vec)
}

#take a df, the figure out {expression, notExpressed}, where notExpress = {low, medium, high} && count(low) == count(medium) == count(high)
assignExprLabels <- function(df,cols){
  labels.df <- as.data.frame(apply(df[as.vector(cols)],2,exprLabels),)
  colnames(labels.df) <-  as.vector(sapply(colnames(labels.df),function(x)paste(x,"ExprLabel",sep=".")))
  cbind(df,labels.df)
}

#for determining the maximal transcript(by sum expression)
rankSumExpr <- function(df){
  df.new <- df[order(df$sumExpr,decreasing=TRUE),]
  df.new$rank <- seq(1,length(df$transcript_id) )
  df.new
}

getExprComp <- function(df){
  with(df,
       data.frame(
         exprInBoth        = length(intersect(which(longPolyAexpr > 0),which(longNonPolyAexpr > 0))),
         none      = length(intersect(which(longPolyAexpr == 0),which(longNonPolyAexpr == 0))),
         either        = length(union(which(longPolyAexpr > 0),which(longNonPolyAexpr > 0))),
         polyA       = length(which(longPolyAexpr > 0)),
         nonPolyA    = length(which(longNonPolyAexpr > 0)),
         polyAOnly   = length(intersect(which(longPolyAexpr > 0),which(longNonPolyAexpr == 0))),
         nonPolyAOnly= length(intersect(which(longPolyAexpr == 0),which(longNonPolyAexpr > 0))),
         total       = length(longNonPolyAexpr)))
}

appendAnalysis <- function(df,cols){
  df <- applyAndAppendFnDf(df,cols,entropy,"entropyExpr")
  df <- applyAndAppendFnDf(df,cols,sum,"sumExpr")
  df <- applyAndAppendFnDf(df,cols,var,"varianceExpr")
  df <- applyAndAppendFnDf(df,cols,mean,"averageExpr")
  df <- applyAndAppendFnDf(df,cols,min,"minExpr")
  df <- applyAndAppendFnDf(df,cols,max,"maxExpr")
  
  df <- applyAndAppendFnDf(df,cols,function(x)maxExpr(x,cols),"maxExprType")
  df <- applyAndAppendFnDf(df=df,cols=cols,FUN=function(x)threeTimesRestSpecifity(x,cols),"threeTimesRest")
  df <- applyAndAppendFnDf(df=df,cols=cols,FUN=function(x)nTimesRestSpecifity(x,cols),"nTimesRest")
  df <- applyAndAppendFnDf(df=df,cols=cols,FUN=function(x)calcTissSpec(x,cols),"tissSpec")
  df <- assignExprLabels(df,cols=cols)
  
  df$expressed = "no"
  df[which(df$maxExprType != "none"),]$expressed <- "yes"
  df
}

getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"transcript_id"]),function(x)x[1]), ~ .id )
applyTranscriptForMaxCols <- function(df,geneCol="gene_id")ldply(split(cd.expr.df,cd.expr.df[c(geneCol)]),getTranscriptForMaxCols)

makeTmpOutFile <- function(x,outdir="/Users/adam/Desktop/"){
  if(missing(x)){
    x = paste(as.vector(proc.time()[1]),"anon",sep="-")
  }
  paste(outdir.tmp,x,sep="")}

#if (!file.exists(outdir)){
#  print("WARNING: cannot access output dir for plots")
#}

getFullGencodeExpr <- function(exprFile = "/Volumes/MyBook/Volume1/scratch/encodeRnaSeq/cshl/allCellsCombined.space",exprCols=2:33){
 appendAnalysis(read.table(file=expr.file,  header=TRUE, stringsAsFactors=FALSE), exprCols)}


preProcessDf<- function(df){
  df <- df[which(df$sumExpr > 0),]
  df$transcript_id <- df$gene_id
  df$Genc_polyA <- rep(0,length(df$gene_id))
  df
}
getPolyACompDf <- function(df, label,exprColIndex){
  df <- preProcessDf(df)
  expr.cols <- sort(colnames(df[exprColIndex]))
  expr.uniq.cols <- unique(unlist(lapply(colnames(df[exprColIndex]),function(x)strsplit(x,".long"))))
  expr.cells = expr.uniq.cols[grep("Poly",expr.uniq.cols,invert=TRUE)]
  expr.cols.polyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longPolyA",sep="."))) 
  expr.cols.nonPolyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longNonPolyA",sep="."))) 
  
 lpa <- melt(df[c("transcript_id",expr.cols.polyA)],measure.vars=sort(expr.cols.polyA),id.vars="transcript_id")
  colnames(lpa) <- c("transcript_id", "expr","longPolyAexpr")
  lpa$cellType <- as.vector(sapply(sapply(as.vector(lpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))
  lpa$seqPullDown <- "longPolyA"
  colnames(lpa) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown" )
  
  lnpa <- melt(df[c("transcript_id",expr.cols.nonPolyA)],measure.vars=sort(expr.cols.nonPolyA),id.vars="transcript_id")
  colnames(lnpa) <- c("transcript_id", "expr","longNonPolyAexpr")
  lnpa$cellType <- as.vector(sapply(sapply(as.vector(lnpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))
  lnpa$seqPullDown <- "longNonPolyA"
  colnames(lnpa) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown")
 
  tm.lnpa <- lnpa[c("transcript_id","cellType","expression")]
  colnames(tm.lnpa) <- c("transcript_id","cellType","longNonPolyA")
  
  tm.lpa <- lpa[c("transcript_id","cellType","expression")]
  colnames(tm.lpa) <- c("transcript_id","cellType","longPolyA")
  
  comb <- merge(tm.lnpa,tm.lpa,by=c("transcript_id","cellType"))
 
  comb <- transform(comb,sum= longNonPolyA +  longPolyA)
  
#   calcMix <- function(lnpa,lpa){
#     sumIn=lnpa+lpa
#     if (sumIn > 0){
#       if (lnpa == 0){
#         "longPolyA"
#       }
#       if (lpa == 0){
#         "longNonPolyA"
#       }
#       else {
#         "both"
#       }
#     }
#       else {
#         "neither"
#       }
#     "default"
#   }
  
  #comb <- transform(comb, expMix = function)
  comb <- transform(comb,logsum= log(longNonPolyA) +  log(longPolyA))  
  comb <- transform(comb,product=longNonPolyA * longPolyA)
  comb$exprMix <- "both"
  comb[intersect(intersect(which(comb$product == 0),which(comb$sum > 0)),which(comb$longNonPolyA == 0)),]$exprMix <- "lpaOnly"
  comb[intersect(intersect(which(comb$product == 0),which(comb$sum > 0)),which(comb$longPolyA == 0)),]$exprMix <- "lnpaOnly"
  comb[intersect(which(comb$product == 0),which(comb$sum == 0)),]$exprMix <- "none"
  mc <- melt(comb,id.var=c("transcript_id","cellType","sum","logsum","exprMix"),measure.va=c("longNonPolyA","longPolyA"))
  colnames(mc)<- c("transcript_id", "cellType", "sum", "logsum", "exprMix","polyApulldown", "RPKM")
  mc$label <- label
  mc
}

compareLncMrnaRnaSeqMix <- function(outdir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/rnaSeq-pullDownCompare",
                                    lnc.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncExprWithStats_transEachSample.tab",
                                    cd.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/cdExprWithStats_transEachSample.tab"){
  
#outdir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/rnaSeq-pullDownCompare"
makeOutFile <- function(x){paste(outdir,x,sep="/")}

#lnc.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncExprWithStats_transEachSample.tab"
lnc.df<-getPolyACompDf(readInTable(lnc.file),"lncRNA",2:33)
lnc.df <- lnc.df[which(lnc.df$sum > 0),]
lnc.df$sumGroups <- cut(lnc.df$sum, breaks=quantile(lnc.df$sum,0:10/10))

#cd.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/cdExprWithStats_transEachSample.tab"
cd.df<-getPolyACompDf(readInTable(cd.file),"mRNA",2:33)
cd.df <- cd.df[which(cd.df$sum > 0),]
cd.df$sumGroups <- cut(cd.df$sum, breaks=quantile(lnc.df$sum,0:10/10))

comb.df <- rbind(cd.df,lnc.df)
comb.df <- transform(comb.df, cellPulldown = paste(cellType,polyApulldown,sep="."))
comb.df <- transform(comb.df, transcriptTypePulldown = paste(label,polyApulldown,sep="."))


makeExprMixForComb <- function(exprMixGroup){
ggplot(comb.df[comb.df$exprMix == exprMixGroup,],aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
  geom_freqpoly(size=1,binwidth=0.1)+
  theme_bw()+
  facet_grid(sumGroups ~ .,scale="free")+
  scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
  ggtitle(paste("transcript pulldown RPKM comparisons\nfacet over bin\nwhere bin=range of transcript expr sum in sample\nTranscripts: ",exprMixGroup,sep=" "))
ggsave( file=makeOutFile(paste("pullDownComp-freqPoly-sumGroups_sampleComb_scaleFree-trans",exprMixGroup,".pdf",sep="")),width=8,height=18)

ggplot(comb.df[comb.df$exprMix == exprMixGroup,],aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
  geom_density()+
  theme_bw()+
  facet_grid(sumGroups ~ .,scale="free")+
  scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
  ggtitle(paste("transcript pulldown RPKM comparisons\nfacet over bin\nwhere bin=range of transcript expr sum in sample\nTranscripts: ",exprMixGroup,sep=" "))
ggsave( file=makeOutFile(paste("pullDownComp-freqPoly-sumGroups_sampleComb_scaleFree-trans",exprMixGroup,"-density.pdf",sep="")),width=8,height=18)

}
sapply(c("both","lnpaOnly","lpaOnly"),makeExprMixForComb)

makeExprMixForCombDistro <- function(exprMixGroup){
  ggplot(comb.df[comb.df$exprMix == exprMixGroup,],aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=1,binwidth=0.1)+
    theme_bw()+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
    ggtitle(paste("transcript pulldown RPKM comparisons\nfacet over bin\nwhere bin=range of transcript expr sum in sample\nTranscripts: ",exprMixGroup,sep=" "))
  ggsave( file=makeOutFile(paste("pullDownComp-freqPoly-totalDistro-trans",exprMixGroup,".pdf",sep="")),width=7,height=7)
 
  ggplot(comb.df[comb.df$exprMix == exprMixGroup,],aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_density()+
    theme_bw()+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
    ggtitle(paste("transcript pulldown RPKM comparisons\nfacet over bin\nwhere bin=range of transcript expr sum in sample\nTranscripts: ",exprMixGroup,sep=" "))
  ggsave( file=makeOutFile(paste("pullDownComp-freqPoly-totalDistro-trans",exprMixGroup,"density.pdf",sep="")),width=7,height=7)
  
}
sapply(c("both","lnpaOnly","lpaOnly"),makeExprMixForCombDistro)

c.df<- comb.df[which(!is.na(comb.df$sumGroups)),]
ggplot(c.df[c.df$exprMix == "both",],aes(x=sumGroups,y=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown))+geom_boxplot(outlier.size=1)+theme_bw()+
  theme( axis.text.x = element_text(angle=25, vjust=0.8))+
  ggtitle("transcriptType-PullDown boxplot for each expression group")+
  xlab("binned expression group")+
  ylab("log(RPKM of transcript)")+
  scale_fill_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+
  scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[5], "mRNA.longPolyA"= brewer.pal(9,"Blues")[9], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[9],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[6]))+
  ggsave( file=makeOutFile("pullDownComp-boxplot-bySumgroups.pdf"),width=10,height=7)


ggplot(c.df, aes(x=log(RPKM),fill=factor(sumGroups)))+geom_histogram()+theme_bw()+
  facet_wrap(polyApulldown~label)+
  ggtitle("Distribution of summation groups in transcript-type/pulldown RPKM distros\ncolor by summation group")+
  xlab("log(RPKM of transcript)")+
  ylab("count")
ggsave( file=makeOutFile("pullDownComp-sumGroupsForPulldown.pdf"),width=10,height=7)

ggplot(c.df, aes(x=log(RPKM),fill=factor(sumGroups)))+geom_histogram(position="fill")+theme_bw()+
  facet_wrap(polyApulldown~label)+
  ggtitle("Distribution of summation groups in transcript-type/pulldown RPKM distros\ncolor by summation group")+
  xlab("log(RPKM of transcript)")+
  ylab("count")
ggsave( file=makeOutFile("pullDownComp-sumGroupsForPulldown-fill.pdf"),width=10,height=7)


rm(c.df)





comb.df
}


compareLncMrnaRnaSeq <- function(outdir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/rnaSeq-pullDownCompare",
                                 lnc.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncExprWithStats_transEachSample.tab",
                                 cd.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/cdExprWithStats_transEachSample.tab"){
  #
  
  #outdir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/rnaSeq-pullDownCompare"
  makeOutFile <- function(x){paste(outdir,x,sep="/")}
  
  #lnc.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncExprWithStats_transEachSample.tab"
  lnc.df<-getPolyACompDf(readInTable(lnc.file),"lncRNA",2:33)
  lnc.df <- lnc.df[which(lnc.df$sum > 0),]
  lnc.df$sumGroups <- cut(lnc.df$sum, breaks=quantile(lnc.df$sum,0:10/10))
  
  #cd.file = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/cdExprWithStats_transEachSample.tab"
  cd.df<-getPolyACompDf(readInTable(cd.file),"mRNA",2:33)
  cd.df <- cd.df[which(cd.df$sum > 0),]
  cd.df$sumGroups <- cut(cd.df$sum, breaks=quantile(lnc.df$sum,0:10/10))
  
  comb.df <- rbind(cd.df,lnc.df)
  ggplot(comb.df,aes(x=log(RPKM),color=polyApulldown,fill=polyApulldown)) + 
    geom_density(alpha=I(0.4))+
    theme_bw()+
    facet_grid(cellType~label)+
    ggtitle("Comparison of mRNA and lncRNA RNA-seq experiments\nFaceted by cell-type and transcript-type")
  ggsave(file=makeOutFile("pullDownComp-density-allBins.pdf"),width=10,height=24)
  
  ggplot(comb.df,aes(x=log(RPKM),color=polyApulldown,fill=polyApulldown)) + 
    geom_freqpoly(size=0.7,binwidth=0.2)+
    theme_bw()+
    facet_grid(cellType~label)+
    ggtitle("Comparison of mRNA and lncRNA RNA-seq experiments\nFaceted by cell-type and transcript-type")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-allBins.pdf"),width=10,height=24)
  
  sapply(levels(cut(lnc.df$sum, breaks=quantile(lnc.df$sum,0:10/10))), 
      function(x.f){
        x.factor <- gsub('\\(','',gsub('\\]','',gsub(",","-",x.f)))
        x.title <- paste("Comparison of mRNA and lncRNA RNA-seq experiments\nFaceted by cell-type and transcript-type\nbin=",x.factor,sep="")
        x.outFile.freq = makeOutFile(paste("pullDownComp-freqPoly-",x.factor,"=binwidth.pdf",sep=""))
        x.outFile.dens = makeOutFile(paste("pullDownComp-freqPoly-",x.factor,"=binwidth.pdf",sep=""))
        ggplot(comb.df,aes(x=log(RPKM),color=polyApulldown,fill=polyApulldown)) + 
          geom_density(alpha=I(0.4))+
          theme_bw()+
          facet_grid(cellType~label)+
          ggtitle(x.title)
        ggsave(file=x.outFile.dens,width=10,height=24)
        
        ggplot(comb.df,aes(x=log(RPKM),color=polyApulldown,fill=polyApulldown)) + 
          geom_freqpoly(size=0.7,binwidth=0.2)+
          theme_bw()+
          facet_grid(cellType~label)+
          ggtitle(x.title)
        ggsave(file=x.outFile.freq,width=10,height=24)
        
      }#end of applied function 
  )#end of sapply over factors...
  
  comb.df <- transform(comb.df, cellPulldown = paste(cellType,polyApulldown,sep="."))
  comb.df <- transform(comb.df, transcriptTypePulldown = paste(label,polyApulldown,sep="."))
  ggplot(comb.df,aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=0.7,binwidth=0.2)+
    theme_bw()+
    facet_grid(cellType~sumGroups,scale="free")+
    scale_fill_manual(values = c("blue","purple","orange","red"))+ 
    scale_color_manual(values = c("blue","purple","orange","red"))+ 
    ggtitle("transcript pulldown RPKM comparisons\nfacet over celltype and bin\nwhere bin=range of transcript expr sum in sample")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-cellTypeVsumGroups.pdf"),width=24,height=24)
  
  
  ggplot(comb.df,aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=0.7,binwidth=0.2)+
    theme_bw()+
    facet_grid(cellType~sumGroups,scale="free_y")+
    scale_fill_manual(values = c(" mRNA.longNonPolyA"="black", "mRNA.longPolyA"="grey", "lncRNA.longPolyA"="blue","lncRNA.longNonPolyA"="red"))+ 

    ggtitle("transcript pulldown RPKM comparisons\nfacet over celltype and bin\nwhere bin=range of transcript expr sum in sample")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-cellTypeVsumGroups_free.pdf"),width=24,height=24)
  
  low.comb.df <- comb.df[which(as.numeric(comb.df$sumGroups) %in% 1:4),]
  ggplot(low.comb.df,aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=1,binwidth=0.3)+
    theme_bw()+
    facet_grid(cellType~sumGroups,scale="free_y")+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
    
    ggtitle("transcript pulldown RPKM comparisons\nfacet over celltype and bin\nwhere bin=range of transcript expr sum in sample")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-cellTypeVsumGroups_low.pdf"),width=24,height=24)
  
  high.comb.df  <- comb.df[which(as.numeric(comb.df$sumGroups) %in% 5:10),]
  ggplot(high.comb.df,aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=0.7,binwidth=0.1)+
    theme_bw()+
    facet_grid(cellType~sumGroups,scale="free_y")+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[7],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[4]))+ 
    
    ggtitle("transcript pulldown RPKM comparisons\nfacet over celltype and bin\nwhere bin=range of transcript expr sum in sample")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-cellTypeVsumGroups_high.pdf"),width=24,height=24)
  
  ggplot(as.data.frame(with(comb.df[which(comb.df$RPKM > 0),],table(sumGroups,transcriptTypePulldown))),aes(x=sumGroups,y=Freq,color=transcriptTypePulldown,fill=transcriptTypePulldown))+
    geom_histogram()+
    theme_bw()+coord_flip()+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[7],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[4]))+ 
    scale_fill_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[3], "mRNA.longPolyA"= brewer.pal(9,"Blues")[7], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[6],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[3]))+ 
    ggtitle("transcript pulldown RPKM counts\nfor each sample by RPKM bin\nOnly transcripts w/ RPKM > 0 counted")
  ggsave(file=makeOutFile("pullDownComp-histogram-forSamples.pdf"),width=7,height=7)
  
  ggplot(comb.df[which(!is.na(comb.df$sumGroups)),],aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=1,binwidth=0.1)+
    theme_bw()+
    facet_grid(sumGroups ~ .,scale="free_y")+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
    
    ggtitle("transcript pulldown RPKM comparisons\nfacet over bin\nwhere bin=range of transcript expr sum in sample")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-sumGroups_sampleComb.pdf"),width=8,height=18)
  
  ggplot(comb.df[which(!is.na(comb.df$sumGroups)),],aes(x=log(RPKM),color=transcriptTypePulldown,fill=transcriptTypePulldown)) + 
    geom_freqpoly(size=1,binwidth=0.1)+
    theme_bw()+
    facet_grid(sumGroups ~ .,scale="free")+
    scale_color_manual(values = c("mRNA.longNonPolyA"= brewer.pal(9,"Blues")[4], "mRNA.longPolyA"= brewer.pal(9,"Blues")[8], "lncRNA.longPolyA"= brewer.pal(9,"Reds")[8],"lncRNA.longNonPolyA"= brewer.pal(9,"Reds")[5]))+ 
    
    ggtitle("transcript pulldown RPKM comparisons\nfacet over bin\nwhere bin=range of transcript expr sum in sample")
  ggsave(file=makeOutFile("pullDownComp-freqPoly-sumGroups_sampleComb_scaleFree.pdf"),width=8,height=18)
  
  
  
  comb.npa.df <- comb.df[which(comb.df$polyApulldown == "longNonPolyA"),]
  comb.npa.df$npa <- comb.npa.df$RPKM
  comb.np.df <- comb.df[which(comb.df$polyApulldown == "longPolyA"),]
 
  comb.npa.df$pa <- comb.pa.df$RPKM
  ggplot(comb.npa.df,aes(x=log(pa),y=log(npa),color=label))+geom_point(size=1)+theme_bw()+facet_grid(cellType ~ sumGroups);
  ggsave(file=makeOutFile("comb-scatterplot.pdf"),height=24,width=24)
  
  
   ggplot(comb.df,aes(x=log(sum),y=log(RPKM),color=polyApulldown))+geom_point(size=1)+facet_wrap(~label,ncol=2)+ theme_bw()
  ggsave(file=makeOutFile("comb-scatterplot-sum-vs-components.pdf"),height=7,width=14)
  
  ggplot(comb.df,aes(x=sum,y=log(RPKM),color=polyApulldown))+stat_smooth(size=1)+facet_wrap(~label,ncol=2)+ theme_bw()
  ggsave(file=makeOutFile("comb-scatterplot-sum-vs-components-qauntile.pdf"),height=7,width=14)
}
  
runAllTransPolyA_Analysis <- function(){
  
  lf <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncExprWithStats_allTrans.tab"
   cd <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/cdExprWithStats_allTrans.tab"
  od <- "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/rnaSeq-pullDownCompare_allTrans"
  compareLncMrnaRnaSeq(outdir=od,lnc.file=lf,cd.file=cd)
  compareLncMrnaRnaSeqMix(outdir=od,lnc.file=lf,cd.file=cd)
}

getAnnotLncDf <- function(lncDf,annotDf,exprCol,annotColName){
  
  lncDf$withinSubset = "false"
  lncDf[which(lncDf$gene_id_short %in% annotDf$gene_id_short),]$withinSubset = "true"
  
  #apply annotDf information to lncDf
  ddply(annotDf,.(gene_id),function(x)x[c("gene_id","lncRnaName")]) -> ensLnc.df
  getLncName <- function(x){ensLnc.df[which(ensLnc.df$gene_id == x ), "lncRnaName"][1]}
  lncDf[c(annotColName)] <- "1"
  lncDf[c(annotColName)] <- apply(lncDf,1,function(x)getLncName(x[c("gene_id_short")]))
  lncDf[is.na(lncDf[c(annotColName)]),annotColName] <- "notFound"
  lncDf
}

makeDir <- function(dir,recursiveCreate=FALSE){
  if (!file.exists(dir)){
    dir.create(path=dir,showWarnings=TRUE,recursive=recursiveCreate,mode="0755")
  }
  
}

getLpaColnames <- function(df,colnamesDf=colnames(df)){
  doubleCols = colnamesDf[as.vector(sapply(colnamesDf,function(x)(typeof(df[1,x]) == "double")))]
  doubleCols[grep("longPolyA$",doubleCols)]
  
}

getLnpaColnames <- function(colnamesDf){
  doubleCols = colnamesDf[as.vector(sapply(colnamesDf,function(x)(typeof(df[1,x]) == "double")))]
  doubleCols[grep("longNonPolyA$",doubleCols)]  
}

editStatsForLncDf <- function(expr.df, cols){
  
  col.names <- colnames(expr.df)[cols]
  expr.df <- applyAndAppendFnDf(expr.df,cols,entropy,"entropyExpr")
  expr.df <- applyAndAppendFnDf(expr.df,cols,sum,"sumExpr")
  expr.df <- applyAndAppendFnDf(expr.df,cols,var,"varianceExpr")
  expr.df <- applyAndAppendFnDf(expr.df,cols,mean,"averageExpr")
  expr.df <- applyAndAppendFnDf(expr.df,cols,min,"minExpr")
  expr.df <- applyAndAppendFnDf(expr.df,cols,max,"maxExpr")
  
  expr.df <- applyAndAppendFnDf(df=expr.df,cols,function(x)maxExpr(x,col.names),"maxExprType")
  expr.df <- applyAndAppendFnDf(df=expr.df,cols=cols,FUN=function(x)threeTimesRestSpecifity(x,col.names),"threeTimesRest")
  expr.df <- applyAndAppendFnDf(df=expr.df,cols=cols,FUN=function(x)nTimesRestSpecifity(x,cols),"nTimesRest")
  expr.df <- applyAndAppendFnDf(df=expr.df,cols=cols,FUN=function(x)calcTissSpec(x,cols),"tissSpec")
  expr.df
  
}

xapplyAndAppendFnDf <-function(df=df,cols=cols,FUN=fn,newCol=newCol){
  df[newCol] <- apply(df[as.vector(cols)],1,function(...){x<-c(...);FUN(x)})
  df
}

xapplyAndAppendFnDfSecondArg <-function(df=df,cols=cols,FUN=fn,newCol=newCol,arg2=arg2){
  df[newCol] <- apply(df[as.vector(cols)],1,function(...){x<-c(...);FUN(x,arg2)})
  df
}

getMemory <- function(){
  gettextf("%.2f Mb stored in memory",
           sum(sapply(unlist(ls(envir=.GlobalEnv)), 
                       function(x)object.size(get(x,envir=.GlobalEnv))))
           / (1000000))
}
getMemory()







