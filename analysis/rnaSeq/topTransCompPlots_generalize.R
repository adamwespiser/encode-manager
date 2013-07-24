home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)


source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R")
#out files created: 
#cd.in.file = getFullPath("data/cdExprWithStats.tab")
#in.file = getFullPath("data/lncExprWithStats.tab")
#combined.in.file = getFullPath("data/combinedExprWithStats.tab")

#outdir<- "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots"
#outdir.tmp <- "/Users/adam/Desktop/top1Trans-tmp-name-";

################### ################### ################### ################### ################### ################### 
################### ################### ################### ################### ################### ################### 


# TEST 
#makeOutFile <- function(x){paste("/Users/adam/test/lncRnaExprT",x,sep="/")}

# IRL
#makeOutFile <- function(x){paste(outdir,x,sep="/")}


################### ################### ################### ################### ################### ################### 
################### ################### ################### ################### ################### ################### 





########################################################
#         LNC RNA 
#         PLOTs
#
########################################################
#

preProcessDf<- function(df){
  df <- df[which(df$sumExpr > 0),]
  df$transcript_id <- df$gene_id
  df$Genc_polyA <- rep(0,length(df$gene_id))
  df
}

runLncTransEachSample <- function(){
combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
lnc.out.dir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots_eachSample"
lnc.fileBase = "lncRNA-gene=MaxTransSample"
lnc.transKeyword = "lncRNA genes as max trans in sample"
lnc.colIndex = 2:33
lnc.df=readInTable(lnc.in.file)

lnc.melt.df <- makePlotsForTranscripts(preProcessDf(lnc.df),
         transKeyword=lnc.transKeyword,
         fileBase=lnc.fileBase,
         outDir=lnc.out.dir,
         exprColIndex=lnc.colIndex,
                    skipPlots=TRUE)  

  
cd.in.file = getFullPath("data/cdExprWithStats_transEachSample.tab")
cd.out.dir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots_eachSample"
cd.fileBase = "mRNA-gene=MaxTransSample"
cd.transKeyword = "mRNA genes as max trans in sample"
cd.colIndex = 2:33
cd.df=readInTable(cd.in.file)
cd.melt.df <-makePlotsForTranscripts(preProcessDf(cd.df),
         transKeyword=cd.transKeyword,
         fileBase=cd.fileBase,
         outDir=cd.out.dir,
                        exprColIndex=cd.colIndex,
                                     skipPlots=TRUE)



makeCombinedPlotsForTranscripts(cd=preProcessDf(cd.df),lnc=preProcessDf(lnc.df)
                        transKeyword="lncRNA & mRNA genes as max trans in sample",
                        fileBase="comb-gene=MaxTransSample",
                        outDir=cd.out.dir,
                        exprColIndex=cd.colIndex,
                          skipPlots=FALSE)



}
#
#
makeCombinedPlotsForTranscripts <- function(cd ,lnc,
                                transKeyword="lncRNA & mRNA genes as max trans in sample",
                                fileBase="comb-gene=MaxTransSample",
                                outDir=cd.out.dir,
                                exprColIndex=cd.colIndex,
                                skipPlots=FALSE){
lnc$transType = "lncRNA"
cd$transType = "mRNA"
if (all.equal(colnames(lnc),colnames(cd))){combined.melt.df <- rbind(lnc,cd)}
printReport <- function(report){ if(verbose == TRUE){print(report)}}
fileBase = paste(outDir,fileBase,sep="/")
makeOutFile <- function(x){paste(fileBase,x,sep="")}

ggplot(combined.melt.df[which(combined.melt.df$variable == "A549.longNonPolyA"),],aes(x=tissSpec,color=transType,fill=transType)) + 
  geom_freqpoly(size=2)+theme_bw()+
  ggtitle("Tissue specificity of lncRNA & mRNA")
ggsave(file=makeOutFile("combined-topTrans-cell-type-tissSpec.pdf"),height=5,width=5)

ggplot(combined.melt.df[which(combined.melt.df$variable == "A549.longNonPolyA"),],aes(x=tissSpec,color=transType,fill=transType)) + 
  geom_density(alpha=I(0.6))+theme_bw()+
  ggtitle("Tissue specificity of lncRNA & mRNA")
ggsave(file=makeOutFile("combined-topTrans-cell-type-tissSpec.pdf.pdf"),height=5,width=5)

ggplot(combined.melt.df,aes(x=log(value),color=transType,fill=transType)) + 
  geom_freqpoly(size=2)+theme_bw()+xlab("log(RPKM)")+
  ggtitle("Total RNA Expression of lncRNA & mRNA\nWhole cell, poly A +/- included")
ggsave(file=makeOutFile("combined-topTrans-cell-type-Expr.pdf"),height=5,width=5)

ggplot(combined.melt.df,aes(x=log(value),color=transType,fill=transType)) + 
  geom_density(alpha=I(0.6))+theme_bw()+xlab("log(RPKM)")+
  ggtitle("Total RNA Expression of lncRNA & mRNA\nWhole cell, poly A +/- included")
ggsave(file=makeOutFile("combined-topTrans-cell-type-Expr-dens.pdf"),height=5,width=5)

ggplot(combined.melt.df,aes(x=log2(value),color=transType))+
  geom_freqpoly(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~variable,ncol=2)+
  ggtitle("cshl transcript lncRNA vs. mRNA expression data\n(Top Expressed Transcript For Each Gene)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("combined-topTrans-cell-type-Expr.pdf"),height=24,width=12)

ggplot(combined.melt.df,aes(x=log2(value),color=transType))+
  geom_freqpoly(binwidth=0.4)+
  theme_bw()+
  ggtitle("cshl transcript lncRNA vs. mRNA expression data\n(Top Expressed Transcript For Each Gene)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("combined-topTrans-Expr.pdf"),height=7,width=7)

as.data.frame(as.list(diag(cor(lnpa.df[expr.cells],lpa.df[expr.cells])))) -> cor
cor$type = "lncRNA"
as.data.frame(as.list(diag(cor(cd.lnpa.df[expr.cells],cd.lpa.df[expr.cells])))) -> cd.cor
cd.cor$type = "mRNA"
comb.cor <- melt(rbind(cd.cor, cor),id.var="type")
colnames(comb.cor) <- c("type","cell","correlationCoeff")
ggplot(comb.cor,aes(x=cell, y=correlationCoeff))+geom_histogram()+facet_wrap(~type,nrow=2)+theme_bw()+coord_flip()+ggtitle("Correlation Coeff between polyA and nonPolyA\nfor each cell type,transcript type")
ggsave(file=makeOutFile("combined-correlation-exprVals.pdf"),height=7,width=7)


cd.ecdf <- ecdf(cd.melt.df$value)
with(cd.df,ecdf(get(expr.cols[1]))) -> col1
qq.col1 <- data.frame(x=as.vector(sapply(cd.df[expr.cols[1]],col1)),y=as.vector(sapply(cd.df[expr.cols[1]],cd.ecdf)))
ggplot(qq.col1,aes(x=x,y=y))+geom_point()+geom_abline(slope=1,intercept=0)+xlim(0,1)+ylim(0,1)
ggsave(makeTmpOutFile("qqplotTest.pdf"),height=7,width=7)

}


#
#


##x in.file
####### TEST EXAMPLE
#df <- readInTable(in.file)
#df <- df[which(df$rank == 1),]
#df <- df[which(df$sumExpr > 0),]
#makePlotsForTranscripts(df=df,transKeyword="lncRNA-testFileOnly",fileBase="fileBaseForFREE",outDir="/Users/adam/test/testDirTop1",exprColIndex=3:34)#

makePlotsForTranscripts <- function(df=df,transKeyword="lncRNA",fileBase="",outDir="/Users/adam/test/",exprColIndex=2:33,verbose=TRUE,skipPlots=FALSE){
  
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  fileBase = paste(outDir,fileBase,sep="/")
  makeOutFile <- function(x){paste(fileBase,x,sep="")}
  
  print("Taking only entries with expression in at least one sample...")
  #df<-df[which(apply(df[,exprColIndex],1,function(...)sum(...)>0)),"sumExpr"]
  
 # if (!any("Genc_polyA" %in% colnames(df))){
#    print("Genc_polyA not detected, setting value to zero")
 #   df$Genc_polyA <- 0}
  
  printReport("Starting data.frame transformation")
  expr.cols <- sort(colnames(df[exprColIndex]))
  expr.uniq.cols <- unique(unlist(lapply(colnames(df[exprColIndex]),function(x)strsplit(x,".long"))))
  expr.cells = expr.uniq.cols[grep("Poly",expr.uniq.cols,invert=TRUE)]
  expr.cols.polyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longPolyA",sep="."))) 
  expr.cols.nonPolyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longNonPolyA",sep="."))) 
  
  lpa <- melt(df[c("transcript_id",expr.cols.polyA)],measure.vars=sort(expr.cols.polyA),id.vars="transcript_id")
  colnames(lpa) <- c("transcript_id", "expr","longPolyAexpr")
  lpa$cellType <- as.vector(sapply(sapply(as.vector(lpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))
  lnpa <- melt(df[c("transcript_id",expr.cols.nonPolyA)],measure.vars=sort(expr.cols.nonPolyA),id.vars="transcript_id")
  colnames(lnpa) <- c("transcript_id", "expr","longNonPolyAexpr")
  lnpa$cellType <- as.vector(sapply(sapply(as.vector(lnpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))
 
  pa.sep.df <- merge(merge(lpa,lnpa,by=c("transcript_id","cellType")),df[c("Genc_polyA","transcript_id")],by="transcript_id")
  lpa$seqPullDown <- "longPolyA"
  colnames(lpa) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown" )
  lnpa$seqPullDown <- "longNonPolyA"
  colnames(lnpa) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown")
  comb <- rbind(lnpa,lpa)
  
  pap.pam.df <- ldply(split(pa.sep.df,pa.sep.df$cellType),getExprComp)
  pam.pam.melt <- melt(pap.pam.df[c(".id","exprInBoth", "none","polyAOnly","nonPolyAOnly")],id.vars=".id")
  pam.pam.melt$variable <- factor(pam.pam.melt$variable,levels=c("exprInBoth","polyAOnly","nonPolyAOnly", "none"))
  lnpa.df<-df[c("transcript_id",expr.cols.nonPolyA)]
  lpa.df<-df[c("transcript_id",expr.cols.polyA)]
  colnames(lpa.df) <- c("transcript_id",expr.cells)
  colnames(lnpa.df) <- c("transcript_id",expr.cells)
  
  printReport("Temp too high...data is melting!!!!")
  melt.df <- melt(df[c("transcript_id",expr.cols,"entropyExpr","threeTimesRest","maxExprType","tissSpec")], 
                       measure.vars=sort(expr.cols),
                       id.vars=c("transcript_id","entropyExpr","threeTimesRest","maxExprType","tissSpec"))
  
  
  melt.df$ctsCat <- "none"
  melt.df[which(melt.df$threeTimesRest == melt.df$variable),]$ctsCat <- "cell-type-this"
  melt.df[intersect(which(melt.df$threeTimesRest != melt.df$variable), which(melt.df$threeTimesRest != "none")),]$ctsCat <- "cell-type-other"
  melt.df$maxExprCat <- "none"
  melt.df[which(melt.df$maxExprType == melt.df$variable),]$maxExprCat <- "cell-type-this"
  melt.df[intersect(which(melt.df$maxExprType != melt.df$variable), which(melt.df$maxExprType != "none")),]$maxExprCat <- "cell-type-other"
  melt.df$jsdCat <- "none"
  melt.df[intersect(which(melt.df$maxExprType == melt.df$variable),which(melt.df$tissSpec > 0.6)),]$jsdCat <- "cell-type-this"
  melt.df[intersect(which(melt.df$maxExprType != melt.df$variable),which(melt.df$tissSpec > 0.6)),]$jsdCat <- "cell-type-other"
  melt.df$tissSpecSample <- "none"
  melt.df$cellTypeSpecificForFacet <- "no"
  melt.df[which(melt.df$variable == melt.df$threeTimesRest),]$cellTypeSpecificForFacet <- "yes"
  melt.df$cellTypeSpecific <- "no"
  melt.df[which(melt.df$threeTimesRest == melt.df$variable),]$cellTypeSpecific <- "yes"
  
  df$isTissSpec <- "notThreeTimes"
  df[which(df$threeTimesRest != "none"),]$isTissSpec <- "threeTimes"
  tmp <- df
  tmp$isTissSpec <- "combinedDistro"
  tissSpec.df <- rbind(df,tmp)
  printReport("Melt Complete...Here by DDPLY")
  mydf_m <- ddply(tissSpec.df, .(isTissSpec), transform, ecd = ecdf(tissSpec)(tissSpec))
  
  maxExprCatCount.df <- ddply(subset(melt.df, value > 0),.(variable,maxExprCat),function(df)dim(df)[1])
  colnames(maxExprCatCount.df) <- c("RnaSeqExpr","ctsCat","transcripts")
  maxExprCatCount.df$ctsCat <- factor(maxExprCatCount.df$ctsCat,levels=c("cell-type-this","cell-type-other","none"))
  
  ctsCount.df <- ddply(subset(melt.df, value > 0),.(variable,ctsCat),function(df)dim(df)[1])
  colnames(ctsCount.df) <- c("RnaSeqExpr","ctsCat","transcripts")
  ctsCount.df$ctsCat <- factor(ctsCount.df$ctsCat,levels=c("cell-type-this","cell-type-other","none"))
  
  jsdCount.df <- ddply(subset(melt.df, value > 0),.(variable,jsdCat),function(df)dim(df)[1])
  colnames(jsdCount.df) <- c("RnaSeqExpr","jsdCat","transcripts")
  jsdCount.df$jsdCat <- factor(jsdCount.df$jsdCat,levels=c("cell-type-this","cell-type-other","none"))
  
  maxExprTable <- as.data.frame(table(df$maxExprType))
  threeTimeTable <- as.data.frame(table(subset(melt.df, threeTimesRest == variable)[c("threeTimesRest")]))
  
  if (skipPlots == TRUE){ return(melt.df) }
  
  printReport("Starting plots...should have just written this in lisp...")
  
  ggplot(pa.sep.df,aes(x=log(longPolyAexpr),y=log(longNonPolyAexpr),color=Genc_polyA))+geom_point(size=1)+
    facet_wrap(~cellType,ncol=3)+
    theme_bw()+
    ggtitle(paste(transKeyword,"\nComparison of cshl longNonPolyA versus longPoly log(RPKM)\nColor is from derrien 2012 data field"))
  ggsave(file=makeOutFile("lncRNA-Compare-PAM-w-PAP.pdf"),height=14,width=7)
  printReport(paste("printed",makeOutFile("lncRNA-Compare-PAM-w-PAP.pdf")))
  ggplot(comb[which(comb$expression > 0),],aes(x=log(expression),fill=seqPullDown))+geom_density(alpha=I(0.4))+
    facet_wrap(~cellType,ncol=3)+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nCell-type expression by RNA-Seq Expr")) +
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorByPullDown.pdf"),height=16,width=8)
  
  ggplot(comb[which(comb$expression > 0),],aes(x=log(expression),color=seqPullDown))+geom_freqpoly(size=1,binwidth=0.5)+
    facet_wrap(~cellType,ncol=3)+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nCell-type expression by RNA-Seq Expr"))+
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorByPullDown-freqpoly.pdf"),height=16,width=8)
  
  ggplot(pam.pam.melt,aes(x=.id,y=value,fill=variable,color=variable))+geom_histogram()+coord_flip()+theme_bw()+
    scale_fill_manual(values = c("exprInBoth" = brewer.pal(9,"Reds")[4],"polyAOnly" = brewer.pal(9,"Greens")[4],"nonPolyAOnly" = brewer.pal(9,"Blues")[4],"none" = brewer.pal(9,"Greys")[2]))+ 
    scale_colour_manual(values = c("exprInBoth" = brewer.pal(9,"Reds")[5],"polyAOnly" = brewer.pal(9,"Greens")[5],"nonPolyAOnly" = brewer.pal(9,"Blues")[5],"none" = brewer.pal(9,"Greys")[3]))+
    ggtitle(paste(transKeyword,"\nTop Transcript Expression:\nComparison of PolyA and NonPolyA pull down experiments\norganized by cell type"))
  ggsave(file=makeOutFile("lncRNA-Compare-PAM-w-PAP-Summary-Chart.pdf"),height=7,width=7)
  
  
  
  
  printReport("Running clustering, call the batman...")
  hclust(dist(as.data.frame(t(as.matrix(subset(df,sumExpr > 0)[expr.cols])))))->hc
  ggdendrogram(hc, rotate=TRUE, size=2)
  ggsave(file=makeOutFile("lncRNA-cellExpr-dendroAll.pdf"),width=12,height=7)
  
  hclust(dist(as.data.frame(t(as.matrix(lnpa.df[expr.cells])))))->hc
  ggdendrogram(hc, rotate=TRUE, size=5)+ggtitle("lncRNA\ntop1 transcript(lncRNA):\npolyA RNA-seq\n")
  ggsave(file=makeOutFile("lncRNA-cellExpr-dendroNonPolyA.pdf"),height=5,width=8)
  
  hclust(dist(as.data.frame(t(as.matrix(lpa.df[expr.cells])))))->hc
  ggdendrogram(hc, rotate=TRUE, size=5)+ggtitle("lncRNA\ntop1 transcript(lncRNA):\nlongPolyA RNA-seq\n")
  ggsave(file=makeOutFile("lncRNA-cellExpr-dendroPolyA.pdf"),height=5,width=8)
  
  ggplot(melt.df[which(melt.df$value>0),],aes(x=log(value),fill=ctsCat))+geom_density(alpha=I(0.4))+
    facet_wrap(~variable,ncol=2,scale="free_y")+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nfilled by category of transcript cell-type specificity of"))+
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorBySpecificity.pdf"),height=24,width=12)
  
  ggplot(melt.df[which(melt.df$value>0),],aes(x=log(value),colour=ctsCat))+geom_freqpoly(size=1,binwidth=0.5)+
    facet_wrap(~variable,ncol=2,scale="free_y")+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nCell-type specificity of transcript is colored")) +
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorBySpecificity-freqPoly.pdf"),height=24,width=12)
  
  ggplot(maxExprCatCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=ctsCat,colour=ctsCat))+geom_histogram()+coord_flip()+ theme_bw()+
    scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
    scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
    ggtitle(paste(transKeyword,"\nCell-type Max Expr\nFor transcripts found in {cell-type/polyA pulldown} pairs"))
  ggsave(file=makeOutFile("cellExpr-maxExprCountHistogram.pdf"),height=7,width=7)
  
  ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=ctsCat,colour=ctsCat))+geom_histogram()+coord_flip()+ theme_bw()+
    scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
    scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
    ggtitle(paste(transKeyword,"\nCell-type Specificity(threeTimesRest)\nFor transcripts found in {cell-type/polyA pulldown} pair"))
  ggsave(file=makeOutFile("cellExpr-specificityCountHistogram.pdf"),height=7,width=7)
  
  ggplot(melt.df[which(melt.df$jsdCat != "none"),],aes(x=variable,y=log2(value),group=transcript_id,color=maxExprType)) + 
    geom_line()+
    theme_bw()+
    facet_wrap(~maxExprType,ncol=2)+
    ggtitle(paste(transKeyword,"\ncshl lncRNA top transcript expression data\nfacet and colory by JSD(0.6) Specificity metric\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("facetJSD.pdf"),height=24,width=16)
  
  ggplot(melt.df[which(melt.df$value>0),],aes(x=log(value),fill=jsdCat))+geom_density(alpha=I(0.4))+
    facet_wrap(~variable,ncol=2,scale="free_y")+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nfilled by category of transcript cell-type JSD(0.6) specificity of"))+
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorByJsdSpecificity.pdf"),height=24,width=12)
  
  ggplot(melt.df[which(melt.df$value>0),],aes(x=log(value),colour=jsdCat))+geom_freqpoly(size=1,binwidth=0.5)+
    facet_wrap(~variable,ncol=2,scale="free_y")+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nCell-type JSD specificity (0.6) of transcript is colored"))+
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorByJSDSpecificity-freqPoly.pdf"),height=24,width=12)
  
  ggplot(melt.df[which(melt.df$value>0),],aes(x=log(value),fill=maxExprCat))+geom_density(alpha=I(0.4))+
    facet_wrap(~variable,ncol=2)+
    ggtitle(paste(transKeyword,"\nCell Type Expression:(top transcripts per gene)\nCell-type max expression of transcript is filled"))+
    theme_bw()
  ggsave(file=makeOutFile("lncRNA-cellExpr-colorByMaxExpr.pdf"),height=24,width=12)
  printReport("just finished 10 plots...")
  #jscCount
  ggplot(jsdCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=jsdCat,colour=jsdCat))+geom_histogram()+coord_flip()+ theme_bw()+
    scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
    scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
    ggtitle(paste(transKeyword,"\nCell-type JSD(0.6) \nFor transcripts found in {cell-type/polyA pulldown} pairs"))
  ggsave(file=makeOutFile("lncRNA-cellExpr-jsdSpecCountHistogram.pdf"),height=7,width=7)
  
  ggplot(maxExprTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
    theme_bw()+
    xlab("Frequency of Cell-type max expression of a lncRNA")+
    ylab("cell-type and extraction method")+
    ggtitle(paste(transKeyword,"\nFrequency of cell types transcript is max expressed in\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("maxExpr-countsPerCell.pdf"),height=7,width=7)
  
  ggplot(melt.df,aes(x=log2(value)))+
    geom_bar(binwidth=0.4)+
    theme_bw()+
    facet_wrap(~variable,ncol=2)+
    ggtitle(paste(transKeyword,"\ncshl transcript lncRNA expression data\n(Top Expressed Transcript For Each Gene)"))+
    xlab("log2(RPKM,whole cell)")+
    ylab("count")
  ggsave(file=makeOutFile("cell-type-Expr.pdf"),height=24,width=16)
  
  threeTimeTable <- as.data.frame(table(df$threeTimesRest))
  ggplot(threeTimeTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
    theme_bw()+
    xlab("Frequency of Cell-type with 3x expression than rest")+
    ylab("cell-type and extraction method")+
    ggtitle(paste(transKeyword,"\nTranscripts in cell-type with 3x RPKM values\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("threeTimeRest-countsPerCell.pdf"),height=7,width=7)
  
  ggplot(threeTimeTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
    theme_bw()+
    xlab("Frequency of Cell-type with 3x expression than rest")+
    ylab("cell-type and extraction method")+
    ggtitle(paste(transKeyword,paste("\nNumber of transcripts in cell-type with 3x RPKM values\n(Top Expressed Transcript For Each Gene)\nN=",
                                     dim(subset(melt.df, threeTimesRest == variable)[c("threeTimesRest")])[1],
                                     "/",
                                     dim(df)[1],
                                     "(total)" )))
  ggsave(file=makeOutFile("threeTimeRest-countsPerCell-noNone.pdf"),height=7,width=7)
  
  ggplot(melt.df,aes(x=log2(value)))+
    geom_bar(binwidth=0.4)+
    theme_bw()+
    facet_wrap(~threeTimesRest,ncol=2)+
    ggtitle(paste(transKeyword,"\ncshl transcript expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types"))+
    xlab("log2(RPKM,whole cell)")+
    ylab("number of transcripts")
  ggsave(file=makeOutFile("threeTimeRest-exprInCellTypes.pdf"),height=24,width=16)
  
  ggplot(melt.df[which(melt.df$threeTimesRest != "none"),],aes(x=log2(value)))+
    geom_bar(binwidth=0.4)+
    theme_bw()+
    facet_wrap(~threeTimesRest,ncol=2)+
    ggtitle(paste(transKeyword,"lncRNA\ncshl transcript lncRNA expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types"))+
    xlab("log2(RPKM,whole cell)")+
    ylab("number of transcripts")
  ggsave(file=makeOutFile("threeTimeRest-exprInCellTypes-NoNone.pdf"),height=24,width=16)
  
  ggplot(melt.df[which(melt.df$threeTimesRest != "none"),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
    geom_line()+
    theme_bw()+
    facet_wrap(~threeTimesRest,ncol=2)+
    ggtitle(paste(transKeyword,paste(transKeyword,"\ncshl  top transcript expression data\nfacet and colory by threeTimesRest Specificity metric\n(Top Expressed Transcript For Each Gene)")))
  ggsave(file=makeOutFile("facetThreeTimesRest.pdf"),height=24,width=16)
  
  ggplot(melt.df[which(melt.df$tissSpec > 0.8),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
    geom_line()+
    theme_bw()+
    facet_wrap(~maxExprType,ncol=2)+
    ggtitle(paste(transKeyword,"\ncshltranscript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.8 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("tissSpecGt0_8-facetMaxExpr.pdf"),height=24,width=16)
  
  ggplot(melt.df[which(melt.df$tissSpec > 0.6),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
    geom_line()+
    theme_bw()+
    facet_wrap(~maxExprType,ncol=2)+
    ggtitle(paste(transKeyword,"\ncshl lncRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("tissSpecGt0_6-facetMaxExpr.pdf"),height=24,width=16)
  printReport("just finished 10 plots...")
  
  ggplot(melt.df[which(melt.df$tissSpec > 0.6),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
    geom_line()+coord_polar()+
    theme_bw()+
    ggtitle(paste(transKeyword,"\ncshl transcript expression data\nfacet by tissue with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("tissSpecGt0_6-combined-polar.pdf"),height=16,width=16)
  
  ggplot(subset(melt.df,threeTimesRest != "none"),aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
    geom_line()+coord_polar()+
    theme_bw()+
    ggtitle(paste(transKeyword,"\ncshl lncRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("threeTimesRest-polar.pdf"),height=16,width=16)
  
  ggplot(melt.df,aes(x=variable,y=log2(value),group=transcript_id,color=cellTypeSpecificForFacet)) + 
    geom_line()+
    theme_bw()+
    facet_wrap(~threeTimesRest,ncol=2)+
    scale_colour_manual(values = c("no"== brewer.pal(9,"Greys")[2],"yes" = brewer.pal(9,"Greys")[8]))+
    ggtitle(paste(transKeyword,"\ncshl transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\n(Top Expressed Transcript For Each Gene)"))
  ggsave(file=makeOutFile("facetThreeTimesRest-colorFacetSpecific.pdf"),height=24,width=16)
  
  ggplot(melt.df[which(melt.df$value > 0),],aes(log(value),color=cellTypeSpecific))+geom_density()+
    theme_bw()+
    xlab("log RPMK")+
    ylab("density")+
    ggtitle(paste(transKeyword,"Cell type-specific expression\nComparing expr in cell type where expression is three time rest\nversus all other instances of  expression\nwithin all cell type"))
  ggsave(file=makeOutFile("lncRNA-cellTypeExpr-Comp-nonCellType.pdf"),height=7,width=7)
  
  #cdf 
  ggplot(mydf_m,aes(x = tissSpec, y = ecd)) + 
    geom_line(aes(group=isTissSpec,colour=isTissSpec,size=3)) +
    theme_bw()+
    ggtitle(paste(transKeyword,"\nCDF of tissSpecifity Measure(JS divergence)\nCummulative Frequency of ThreeTimesRest \nvs. Non-ThreeTimes rest detected"))
  ggsave(file=makeOutFile("lncRNA-tissSpec-CDF.pdf"),height=7,width=7)
  
  
  if(!("nTimesRest" %in% colnames(df))){
    df <- applyAndAppendFnDf(df=df,cols=expr.cols,FUN=function(x)nTimesRestSpecifity(x,expr.cols),"nTimesRest")
  }
  
  ggplot(df, aes(x=(tissSpec),y=log(nTimesRest)))+geom_point(alpha=I(0.4),size=1)+
    theme_bw()+
    geom_line(aes(y=log(3),color="3 times rest"))+
    geom_line(aes(y=log(5),color="5 times rest"))+
    layer(
      data=data.frame(x=seq(0,1,0.001)),
      mapping=aes(x=x,y=x),
      color="green",
      geom="line",
      size=0
    )+
    #geom_abline(intercept = 0, slope = 1)
    ggtitle(paste(transKeyword,paste("Comparison of tissue specificity \nand (top expression n times greater than next)\nN=",dim(df)[1])))
  ggsave(file=makeOutFile("lncRNA-tissSpec-vs-nTimesRest.pdf"),height=7,width=7)
  printReport("Function Complete. Pass go and collect $200.")

  
  melt.df
}
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
##         END OF FUNCTION
###
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

###
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

###
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

###
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################

###
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################
#############################################################################################################################################################################################




##############################################################################################################################################################################################
############################################################################################################################################################################################
#
#    combined lncRNA/mRNA
#
#
##############################################################################################################################################################################################################

###########################
#     CELL/PULLDOWN
#     EXPR 
########################### ggplot(comb[which(comb$expression > 0),],aes(x=log(expression),color=seqPullDown))+geom_freqpoly(size=1,binwidth=0.5)+


