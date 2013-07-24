home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)
library(hexbin)
library(RColorBrewer)
library(reshape)
library(ggplot2)
source(getFullPath("analysis/RnaSeq/exprLib.R"))

#out files created: 
cd.in.file = getFullPath("data/cdExprWithStats.tab")
lnc.in.file = getFullPath("data/lncExprWithStats.tab")
combined.in.file = getFullPath("data/combinedExprWithStats.tab")

outdir<- "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots"
outdir.tmp <- "/Users/adam/Desktop/top1Trans-tmp-name-";

################### ################### ################### ################### ################### ################### 
################### ################### ################### ################### ################### ################### 


# TEST 
makeOutFile <- function(x){paste("/Users/adam/test/lncRnaExprT",x,sep="/")}

# IRL
#makeOutFile <- function(x){paste(outdir,x,sep="/")}


################### ################### ################### ################### ################### ################### 
################### ################### ################### ################### ################### ################### 



makeTmpOutFile <- function(x){
              if(missing(x)){
                x = paste(as.vector(proc.time()[1]),"anon",sep="-")
              }
              paste(outdir.tmp,x,sep="")}

if (!file.exists(outdir)){
  print("WARNING: cannot access output dir for plots")
}


########################################################
#         LNC RNA 
#         PLOTs
#
########################################################
lnc.df <- readInTable(lnc.in.file)
expr.cols <- colnames(lnc.df[,3:34])
expr.uniq.cols <- unique(unlist(lapply(colnames(lnc.df[,3:34]),function(x)strsplit(x,".long"))))
expr.cells = expr.uniq.cols[grep("Poly",expr.uniq.cols,invert=TRUE)]
expr.cols.polyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longPolyA",sep="."))) 
expr.cols.nonPolyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longNonPolyA",sep="."))) 
lnc.top1.df <- lnc.df[which(lnc.df$rank == 1),]
lnc.top1.df <- lnc.top1.df[which(lnc.top1.df$sumExpr > 0),]

lnc.melt.top1.df <- melt(lnc.top1.df[c("transcript_id",expr.cols,"entropyExpr","threeTimesRest","maxExprType","tissSpec")], 
                         measure.vars=sort(expr.cols),
                         id.vars=c("transcript_id","entropyExpr","threeTimesRest","maxExprType","tissSpec"))


lnc.lpa <- melt(lnc.top1.df[c("transcript_id",expr.cols.polyA)],measure.vars=sort(expr.cols.polyA),id.vars="transcript_id")
colnames(lnc.lpa) <- c("transcript_id", "expr","longPolyAexpr")
lnc.lpa$cellType <- as.vector(sapply(sapply(as.vector(lnc.lpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))

lnc.lnpa <- melt(lnc.top1.df[c("transcript_id",expr.cols.nonPolyA)],measure.vars=sort(expr.cols.nonPolyA),id.vars="transcript_id")
colnames(lnc.lnpa) <- c("transcript_id", "expr","longNonPolyAexpr")
lnc.lnpa$cellType <- as.vector(sapply(sapply(as.vector(lnc.lnpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))

lnc.pa.sep.df <- merge(lnc.lpa,lnc.lnpa,by=c("transcript_id","cellType"))
lnc.pa.sep.df <- merge(lnc.pa.sep.df,lnc.top1.df[c("Genc_polyA","transcript_id")],by="transcript_id")
ggplot(lnc.pa.sep.df,aes(x=log(longPolyAexpr),y=log(longNonPolyAexpr),color=Genc_polyA))+geom_point(size=1)+
  facet_wrap(~cellType,ncol=3)+
  theme_bw()+
  ggtitle("lncRNA\nComparison of cshl longNonPolyA versus longPoly log(RPKM)\nColor is from derrien 2012 data field")
ggsave(file=makeOutFile("lncRNA-Compare-PAM-w-PAP.pdf"),height=14,width=7)

lnc.lpa$seqPullDown <- "longPolyA"
lp.df <- lnc.lpa
colnames(lp.df) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown" )
lnc.lnpa$seqPullDown <- "longNonPolyA"
lnp.df <- lnc.lnpa
colnames(lnp.df) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown")
lnc.comb <- rbind(lnp.df,lp.df)

ggplot(lnc.comb[which(lnc.comb$expression > 0),],aes(x=log(expression),fill=seqPullDown))+geom_density(alpha=I(0.4))+
  facet_wrap(~cellType,ncol=3)+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nCell-type expression by RNA-Seq Expr")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorByPullDown.pdf"),height=16,width=8)

ggplot(lnc.comb[which(lnc.comb$expression > 0),],aes(x=log(expression),color=seqPullDown))+geom_freqpoly(size=1,binwidth=0.5)+
  facet_wrap(~cellType,ncol=3)+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nCell-type expression by RNA-Seq Expr")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorByPullDown-freqpoly.pdf"),height=16,width=8)



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
lnc.pap.pam.df <- ldply(split(lnc.pa.sep.df,lnc.pa.sep.df$cellType),getExprComp)
lnc.pam.pam.melt <- melt(lnc.pap.pam.df[c(".id","exprInBoth", "none","polyAOnly","nonPolyAOnly")],id.vars=".id")


lnc.pam.pam.melt$variable <- factor(lnc.pam.pam.melt$variable,levels=c("exprInBoth","polyAOnly","nonPolyAOnly", "none"))
ggplot(lnc.pam.pam.melt,aes(x=.id,y=value,fill=variable,color=variable))+geom_histogram()+coord_flip()+theme_bw()+
  scale_fill_manual(values = c("exprInBoth" = brewer.pal(9,"Reds")[4],"polyAOnly" = brewer.pal(9,"Greens")[4],"nonPolyAOnly" = brewer.pal(9,"Blues")[4],"none" = brewer.pal(9,"Greys")[2]))+ 
  scale_colour_manual(values = c("exprInBoth" = brewer.pal(9,"Reds")[5],"polyAOnly" = brewer.pal(9,"Greens")[5],"nonPolyAOnly" = brewer.pal(9,"Blues")[5],"none" = brewer.pal(9,"Greys")[3]))+
  ggtitle("lncRNA\nTop Transcript Expression:\nComparison of PolyA and NonPolyA pull down experiments\norganized by cell type")
ggsave(file=makeOutFile("lncRNA-Compare-PAM-w-PAP-Summary-Chart.pdf"),height=7,width=7)



#dendro
lnc.lnpa.df<-lnc.top1.df[c("transcript_id",expr.cols.nonPolyA)]
lnc.lpa.df<-lnc.top1.df[c("transcript_id",expr.cols.polyA)]
colnames(lnc.lpa.df) <- c("transcript_id",expr.cells)
colnames(lnc.lnpa.df) <- c("transcript_id",expr.cells)


hclust(dist(as.data.frame(t(as.matrix(subset(lnc.top1.df,sumExpr > 0)[expr.cols])))))->lnc.hc
ggdendrogram(lnc.hc, rotate=TRUE, size=2)
ggsave(file=makeOutFile("lncRNA-cellExpr-dendroAll.pdf"),width=12)
 
 hclust(dist(as.data.frame(t(as.matrix(lnc.lnpa.df[expr.cells])))))->lnc.hc
ggdendrogram(lnc.hc, rotate=TRUE, size=5)+ggtitle("lncRNA\ntop1 transcript(lncRNA):\npolyA RNA-seq\n")
 ggsave(file=makeOutFile("lncRNA-cellExpr-dendroNonPolyA.pdf"),height=5,width=8)
 
hclust(dist(as.data.frame(t(as.matrix(lnc.lpa.df[expr.cells])))))->lnc.hc
ggdendrogram(lnc.hc, rotate=TRUE, size=5)+ggtitle("lncRNA\ntop1 transcript(lncRNA):\nlongPolyA RNA-seq\n")
ggsave(file=makeOutFile("lncRNA-cellExpr-dendroPolyA.pdf"),height=5,width=8)


lnc.melt.top1.df$ctsCat <- "none"
lnc.melt.top1.df[which(lnc.melt.top1.df$threeTimesRest == lnc.melt.top1.df$variable),]$ctsCat <- "cell-type-this"
lnc.melt.top1.df[intersect(which(lnc.melt.top1.df$threeTimesRest != lnc.melt.top1.df$variable), which(lnc.melt.top1.df$threeTimesRest != "none")),]$ctsCat <- "cell-type-other"

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$value>0),],aes(x=log(value),fill=ctsCat))+geom_density(alpha=I(0.4))+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nfilled by category of transcript cell-type specificity of")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorBySpecificity.pdf"),height=24,width=12)


ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$value>0),],aes(x=log(value),colour=ctsCat))+geom_freqpoly(size=1,binwidth=0.5)+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nCell-type specificity of transcript is colored")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorBySpecificity-freqPoly.pdf"),height=24,width=12)

lnc.melt.top1.df$maxExprCat <- "none"
lnc.melt.top1.df[which(lnc.melt.top1.df$maxExprType == lnc.melt.top1.df$variable),]$maxExprCat <- "cell-type-this"
lnc.melt.top1.df[intersect(which(lnc.melt.top1.df$maxExprType != lnc.melt.top1.df$variable), which(lnc.melt.top1.df$maxExprType != "none")),]$maxExprCat <- "cell-type-other"

ctsCount.df <- ddply(subset(lnc.melt.top1.df, value > 0),.(variable,maxExprCat),function(df)dim(df)[1])
colnames(ctsCount.df) <- c("RnaSeqExpr","ctsCat","transcripts")
ctsCount.df$ctsCat <- factor(ctsCount.df$ctsCat,levels=c("cell-type-this","cell-type-other","none"))
ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=ctsCat,colour=ctsCat))+geom_histogram()+coord_flip()+ theme_bw()+
  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
  ggtitle("lncRNA\nCell-type Max Expr\nFor transcripts found in {cell-type/polyA pulldown} pairs")
ggsave(file=makeOutFile("lncRNA-cellExpr-maxExprCountHistogram.pdf"),height=7,width=7)

ctsCount.df <- ddply(subset(lnc.melt.top1.df, value > 0),.(variable,ctsCat),function(df)dim(df)[1])
colnames(ctsCount.df) <- c("RnaSeqExpr","ctsCat","transcripts")
ctsCount.df$ctsCat <- factor(ctsCount.df$ctsCat,levels=c("cell-type-this","cell-type-other","none"))
ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=ctsCat,colour=ctsCat))+geom_histogram()+coord_flip()+ theme_bw()+
  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
  ggtitle("lncRNA\nCell-type Specificity(threeTimesRest)\nFor transcripts found in {cell-type/polyA pulldown} pair")
ggsave(file=makeOutFile("lncRNA-cellExpr-specificityCountHistogram.pdf"),height=7,width=7)


####### JSD
lnc.melt.top1.df$jsdCat <- "none"
lnc.melt.top1.df[intersect(which(lnc.melt.top1.df$maxExprType == lnc.melt.top1.df$variable),which(lnc.melt.top1.df$tissSpec > 0.6)),]$jsdCat <- "cell-type-this"
lnc.melt.top1.df[intersect(which(lnc.melt.top1.df$maxExprType != lnc.melt.top1.df$variable),which(lnc.melt.top1.df$tissSpec > 0.6)),]$jsdCat <- "cell-type-other"
lnc.melt.top1.df$tissSpecSample <- "none"

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$jsdCat != "none"),],aes(x=variable,y=log2(value),group=transcript_id,color=maxExprType)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("lncRNA\ncshl lncRNA top transcript expression data\nfacet and colory by JSD(0.6) Specificity metric\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-facetJSD.pdf"),height=24,width=16)





ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$value>0),],aes(x=log(value),fill=jsdCat))+geom_density(alpha=I(0.4))+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nfilled by category of transcript cell-type JSD(0.6) specificity of")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorByJsdSpecificity.pdf"),height=24,width=12)


ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$value>0),],aes(x=log(value),colour=jsdCat))+geom_freqpoly(size=1,binwidth=0.5)+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nCell-type JSD specificity (0.6) of transcript is colored")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorByJSDSpecificity-freqPoly.pdf"),height=24,width=12)


jsdCount.df <- ddply(subset(lnc.melt.top1.df, value > 0),.(variable,jsdCat),function(df)dim(df)[1])
colnames(jsdCount.df) <- c("RnaSeqExpr","jsdCat","transcripts")
jsdCount.df$jsdCat <- factor(jsdCount.df$jsdCat,levels=c("cell-type-this","cell-type-other","none"))
ggplot(jsdCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=jsdCat,colour=jsdCat))+geom_histogram()+coord_flip()+ theme_bw()+
  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
  ggtitle("lncRNA\nCell-type JSD(0.6) \nFor transcripts found in {cell-type/polyA pulldown} pairs")
ggsave(file=makeOutFile("lncRNA-cellExpr-jsdSpecCountHistogram.pdf"),height=7,width=7)

###END JS

#######begin edit
#cd.melt.top1.df$maxExprCat <- "none"
#cd.melt.top1.df[which(cd.melt.top1.df$maxExprType == cd.melt.top1.df$variable),]$maxExprCat <- "cell-type-this"
#cd.melt.top1.df[intersect(which(cd.melt.top1.df$maxExprType != cd.melt.top1.df$variable), which(cd.melt.top1.df$maxExprType != "none")),]$maxExprCat <- "cell-type-other"

#ctsCount.df <- ddply(subset(cd.melt.top1.df, value > 0),.(variable,maxExprCat),function(df)dim(df)[1])
#colnames(ctsCount.df) <- c("RnaSeqExpr","maxExprCat","transcripts")
#ctsCount.df$maxExprCat <- factor(ctsCount.df$maxExprCat,levels=c("cell-type-this","cell-type-other","none"))
#ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=maxExprCat,colour=maxExprCat))+geom_histogram()+coord_flip()+ theme_bw()+
#  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
#  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
#  ggtitle("Cell-type Max Expressed Transcript\nFor mRNA transcripts found in {cell-type/polyA pulldown} pair")
#ggsave(file=makeOutFile("mRNA-cellExpr-maxExprCountHistogram.pdf"),height=7,width=7)

#cd.melt.top1.df$maxSpecCat <- "none"
#cd.melt.top1.df[which(cd.melt.top1.df$threeTimesRest == cd.melt.top1.df$variable),]$maxSpecCat <- "cell-type-this"
#cd.melt.top1.df[intersect(which(cd.melt.top1.df$threeTimesRest != cd.melt.top1.df$variable), which(cd.melt.top1.df$threeTimesRest != "none")),]$maxSpecCat <- "cell-type-other"

#ctsCount.df <- ddply(subset(cd.melt.top1.df, value > 0),.(variable,maxSpecCat),function(df)dim(df)[1])
#colnames(ctsCount.df) <- c("RnaSeqExpr","maxSpecCat","transcripts")
#ctsCount.df$maxExprCat <- factor(ctsCount.df$maxSpecCat,levels=c("cell-type-this","cell-type-other","none"))
#ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=maxSpecCat,colour=maxSpecCat))+geom_histogram()+coord_flip()+ theme_bw()+
#  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
#  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
#  ggtitle("Cell-type Specific(ThreeTimesRest) Expressed Transcript\nFor mRNA transcripts found in {cell-type/polyA pulldown} pair")

#ggsave(file=makeOutFile("mRNA-cellExpr-specifityCountHistogram.pdf"),height=7,width=7)
###########end edit


ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$value>0),],aes(x=log(value),fill=maxExprCat))+geom_density(alpha=I(0.4))+
  facet_wrap(~variable,ncol=2)+
  ggtitle("lncRNA\nCell Type Expression:(top transcripts per gene)\nCell-type max expression of transcript is filled")+
  theme_bw()
ggsave(file=makeOutFile("lncRNA-cellExpr-colorByMaxExpr.pdf"),height=24,width=12)

maxExprTable <- as.data.frame(table(lnc.top1.df$maxExprType))
ggplot(maxExprTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
  theme_bw()+
  xlab("Frequency of Cell-type max expression of a lncRNA")+
  ylab("cell-type and extraction method")+
  ggtitle("lncRNA\nFrequency of cell types transcript is max expressed in\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-maxExpr-countsPerCell.pdf"),height=7,width=7)


ggplot(lnc.melt.top1.df,aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~variable,ncol=2)+
  ggtitle("lncRNA\ncshl transcript lncRNA expression data\n(Top Expressed Transcript For Each Gene)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("lncRNA-topTrans-cell-type-Expr.pdf"),height=24,width=16)

threeTimeTable <- as.data.frame(table(lnc.top1.df$threeTimesRest))
ggplot(threeTimeTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
  theme_bw()+
  xlab("Frequency of Cell-type with 3x expression than rest")+
  ylab("cell-type and extraction method")+
  ggtitle("lncRNA\nTranscripts in cell-type with 3x RPKM values\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-threeTimeRest-countsPerCell.pdf"))


threeTimeTable <- as.data.frame(table(subset(lnc.melt.top1.df, threeTimesRest == variable)[c("threeTimesRest")]))
ggplot(threeTimeTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
  theme_bw()+
  xlab("Frequency of Cell-type with 3x expression than rest")+
  ylab("cell-type and extraction method")+
  ggtitle(paste("lncRNA\nNumber of transcripts in cell-type with 3x RPKM values\n(Top Expressed Transcript For Each Gene)\nN=",
                dim(subset(lnc.melt.top1.df, threeTimesRest == variable)[c("threeTimesRest")])[1],
                "/",
                dim(lnc.top1.df)[1],
                "(total)" ))
ggsave(file=makeOutFile("lncRNA-topTrans-threeTimeRest-countsPerCell-noNone.pdf"))



ggplot(lnc.melt.top1.df,aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  ggtitle("lncRNA\ncshl transcript lncRNA expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types")+
  xlab("log2(RPKM,whole cell)")+
  ylab("number of transcripts")
ggsave(file=makeOutFile("lncRNA-topTrans-threeTimeRest-exprInCellTypes.pdf"),height=24,width=16)

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$threeTimesRest != "none"),],aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  ggtitle("lncRNA\ncshl transcript lncRNA expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types")+
  xlab("log2(RPKM,whole cell)")+
  ylab("number of transcripts")
ggsave(file=makeOutFile("lncRNA-topTrans-threeTimeRest-exprInCellTypes-NoNone.pdf"),height=24,width=16)

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$threeTimesRest != "none"),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  ggtitle("lncRNA\ncshl lncRNA top transcript expression data\nfacet and colory by threeTimesRest Specificity metric\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-facetThreeTimesRest.pdf"),height=24,width=16)

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$tissSpec > 0.8),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("lncRNA\ncshltranscript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.8 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-tissSpecGt0_8-facetMaxExpr.pdf"),height=24,width=16)

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$tissSpec > 0.6),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("lncRNA\ncshl lncRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-tissSpecGt0_6-facetMaxExpr.pdf"),height=24,width=16)


#polar
ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$tissSpec > 0.6),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+coord_polar()+
  theme_bw()+
  ggtitle("lncRNA\ncshl transcript expression data\nfacet by tissue with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-tissSpecGt0_6-combined-polar.pdf"),height=16,width=16)

# TOO LONG TO RUN
#ggplot(lnc.melt.top1.df,aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
#  geom_line()+coord_polar()+
#  theme_bw()+
#  ggtitle("cshl lncRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity\n(Top Expressed Transcript For Each Gene)")
#ggsave(file=makeOutFile("lncRNA-topTrans-combined-polar.pdf"),height=16,width=16)

ggplot(subset(lnc.melt.top1.df,threeTimesRest != "none"),aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+coord_polar()+
  theme_bw()+
  ggtitle("lncRNA\ncshl lncRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity(max JSD)\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-threeTimesRest-polar.pdf"),height=16,width=16)


lnc.melt.top1.df$cellTypeSpecificForFacet <- "no"
lnc.melt.top1.df[which(lnc.melt.top1.df$variable == lnc.melt.top1.df$threeTimesRest),]$cellTypeSpecificForFacet <- "yes"
ggplot(lnc.melt.top1.df,aes(x=variable,y=log2(value),group=transcript_id,color=cellTypeSpecificForFacet)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  scale_colour_manual(values = c("no"== brewer.pal(9,"Greys")[2],"yes" = brewer.pal(9,"Greys")[8]))+
  ggtitle("lncRNA\ncshl transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("lncRNA-topTrans-facetThreeTimesRest-colorFacetSpecific.pdf"),height=24,width=16)




#princomp(lnc.top1.df[expr.cols])->lnc.pca
#as.data.frame(lnc.pca$scores) -> lnc.pca.factors
#lnc.pca.factors$transcript_id <- lnc.top1.df$transcript_id
#lnc.pca.melt <- melt(lnc.pca.factors)
#ggplot(lnc.pca.factors,aes(x=Comp.1,fill=Comp.2))+geom_density()
#ggsave("~/Desktop/pcaCompsForLncTopGene.pdf")


lnc.melt.top1.df$cellTypeSpecific <- "no"
lnc.melt.top1.df[which(lnc.melt.top1.df$threeTimesRest == lnc.melt.top1.df$variable),]$cellTypeSpecific <- "yes"

ggplot(lnc.melt.top1.df[which(lnc.melt.top1.df$value > 0),],aes(log(value),color=cellTypeSpecific))+geom_density()+
  theme_bw()+
  xlab("log RPMK")+
  ylab("density")+
  ggtitle("Cell type-specific expression\nComparing lncRNA in cell type where expression is three time rest\nversus all other instances of lncRNA expression\nwithin all cell type")
ggsave(file=makeOutFile("lncRNA-cellTypeExpr-Comp-nonCellType.pdf"))


#cdf 
lnc.top1.df$isTissSpec <- "notThreeTimes"
lnc.top1.df[which(lnc.top1.df$threeTimesRest != "none"),]$isTissSpec <- "threeTimes"
tmp <- lnc.top1.df
tmp$isTissSpec <- "combinedDistro"
tissSpec.df <- rbind(lnc.top1.df,tmp)
#subset(lnc.top1.df, threeTimesRest != "none") -> lnc.spec.df
mydf_m <- ddply(tissSpec.df, .(isTissSpec), transform, ecd = ecdf(tissSpec)(tissSpec))
#mydf_m <- ddply(lnc.spec.df, .(variable), transform, ecd = ecdf(tissSpec)(tissSpec))

ggplot(mydf_m,aes(x = tissSpec, y = ecd)) + 
  geom_line(aes(group=isTissSpec,colour=isTissSpec,size=3)) +
  theme_bw()+
  ggtitle("CDF of tissSpecifity Measure(JS divergence)\nCummulative Frequency of ThreeTimesRest \nvs. Non-ThreeTimes rest detected(lncRNA)")
ggsave(file=makeOutFile("lncRNA-tissSpec-CDF.pdf"),height=7,width=7)


if(!("nTimesRest" %in% colnames(lnc.top1.df))){
lnc.top1.df <- applyAndAppendFnDf(df=lnc.top1.df,cols=expr.cols,FUN=function(x)nTimesRestSpecifity(x,expr.cols),"nTimesRest")
}

ggplot(lnc.top1.df, aes(x=(tissSpec),y=log(nTimesRest)))+geom_point(alpha=I(0.4),size=1)+
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
  ggtitle(paste("Comparison of tissue specificity \nand (top expression n times greater than next)\nN=",dim(lnc.top1.df)[1]))
ggsave(file=makeOutFile("lncRNA-tissSpec-vs-nTimesRest.pdf"),height=7,width=7)
########################################################
#         messenger RNA 
#         PLOTs
#
########################################################

cd.df <- readInTable(cd.in.file)
cd.top1.df <- cd.df[which(cd.df$rank == 1),]
cd.top1.df <- cd.top1.df[which(cd.top1.df$sumExpr > 0),]
expr.cols <- colnames(lnc.df[,3:34])
rm(list=c(cd.df,lnc.df))
cd.melt.top1.df <- melt(cd.top1.df[c("transcript_id",expr.cols,"entropyExpr","threeTimesRest","maxExprType","tissSpec")], 
                         measure.vars=sort(expr.cols),
                         id.vars=c("transcript_id","entropyExpr","threeTimesRest","maxExprType","tissSpec"))


pc.lpa <- melt(cd.top1.df[c("transcript_id",expr.cols.polyA)],measure.vars=sort(expr.cols.polyA),id.vars="transcript_id")
colnames(pc.lpa) <- c("transcript_id", "expr","longPolyAexpr")
pc.lpa$cellType <- as.vector(sapply(sapply(as.vector(pc.lpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))

pc.lnpa <- melt(cd.top1.df[c("transcript_id",expr.cols.nonPolyA)],measure.vars=sort(expr.cols.nonPolyA),id.vars="transcript_id")
colnames(pc.lnpa) <- c("transcript_id", "expr","longNonPolyAexpr")
pc.lnpa$cellType <- as.vector(sapply(sapply(as.vector(pc.lnpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))

pc.pa.sep.df <- merge(pc.lpa,pc.lnpa,by=c("transcript_id","cellType"))
#pc.pa.sep.df <- merge(pc.pa.sep.df,pc.top1.df[c("Genc_polyA","transcript_id")],by="transcript_id")
ggplot(pc.pa.sep.df,aes(x=log(longPolyAexpr),y=log(longNonPolyAexpr)))+geom_point(size=1)+
  facet_wrap(~cellType,ncol=3)+
  theme_bw()+
  ggtitle("Comparison of cshl longNonPolyA versus longPoly log(RPKM)")
ggsave(file=makeOutFile("mRNA-Compare-PAM-w-PAP.pdf"),height=14,width=7)

pc.lpa$seqPullDown <- "longPolyA"
pc.lp.df <- pc.lpa
colnames(pc.lp.df) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown" )
pc.lnpa$seqPullDown <- "longNonPolyA"
pc.lnp.df <- pc.lnpa
colnames(pc.lnp.df) <- c("transcript_id", "expr",  "expression" ,"cellType" ,"seqPullDown")
pc.comb <- rbind(pc.lnp.df,pc.lp.df)

ggplot(pc.comb[which(pc.comb$expression > 0),],aes(x=log(expression),fill=seqPullDown))+geom_density(alpha=I(0.4))+
  facet_wrap(~cellType,ncol=3)+
  ggtitle("Cell Type Expression:(top mRNA transcripts per gene)\nCell-type expression by RNA-Seq Expr")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorByPullDown.pdf"),height=16,width=8)

ggplot(pc.comb[which(pc.comb$expression > 0),],aes(x=log(expression),color=seqPullDown))+geom_freqpoly(size=1,binwidth=0.5)+
  facet_wrap(~cellType,ncol=3)+
  ggtitle("Cell Type Expression:(top mRNA transcripts per gene)\nCell-type expression by RNA-Seq Expr")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorByPullDown-freqpoly.pdf"),height=16,width=8)




pc.pap.pam.df <- ldply(split(pc.pa.sep.df,pc.pa.sep.df$cellType),getExprComp)
pc.pam.pam.melt <- melt(pc.pap.pam.df[c(".id","exprInBoth", "none","polyAOnly","nonPolyAOnly")],id.vars=".id")
pc.pam.pam.melt$variable <- factor(pc.pam.pam.melt$variable,levels=c("exprInBoth","polyAOnly","nonPolyAOnly", "none"))

ggplot(pc.pam.pam.melt,aes(x=.id,y=value,fill=variable,color=variable))+geom_histogram()+coord_flip()+theme_bw()+
  scale_fill_manual(values = c("exprInBoth" = brewer.pal(9,"Reds")[4],"polyAOnly" = brewer.pal(9,"Greens")[4],"nonPolyAOnly" = brewer.pal(9,"Blues")[4],"none" = brewer.pal(9,"Greys")[2]))+ 
  scale_colour_manual(values = c("exprInBoth" = brewer.pal(9,"Reds")[5],"polyAOnly" = brewer.pal(9,"Greens")[5],"nonPolyAOnly" = brewer.pal(9,"Blues")[5],"none" = brewer.pal(9,"Greys")[3]))+
  ggtitle("Top Transcript Expression:\nComparison of PolyA and NonPolyA pull down experiments")
ggsave(file=makeOutFile("mRNA-Compare-PAM-w-PAP-Summary-Chart.pdf"),height=7,width=7)



#dendro
cd.lnpa.df<-cd.top1.df[c("transcript_id",expr.cols.nonPolyA)]
cd.lpa.df<-cd.top1.df[c("transcript_id",expr.cols.polyA)]
colnames(cd.lpa.df) <- c("transcript_id",expr.cells)
colnames(cd.lnpa.df) <- c("transcript_id",expr.cells)


hclust(dist(as.data.frame(t(as.matrix(subset(cd.top1.df,sumExpr > 0)[expr.cols])))))->cd.hc
ggdendrogram(cd.hc, rotate=TRUE, size=2)
ggsave(file=makeOutFile("mRNA-cellExpr-dendroAll.pdf"),width=12)

hclust(dist(as.data.frame(t(as.matrix(cd.lnpa.df[expr.cells])))))->cd.hc
ggdendrogram(cd.hc, rotate=TRUE, size=5)+ggtitle("top1 transcript(mRNA):\npolyA RNA-seq\n")
ggsave(file=makeOutFile("mRNA-cellExpr-dendroNonPolyA.pdf"),height=5,width=8)

hclust(dist(as.data.frame(t(as.matrix(cd.lpa.df[expr.cells])))))->cd.hc
ggdendrogram(cd.hc, rotate=TRUE, size=5)+ggtitle("top1 transcript(mRNA):\nlongPolyA RNA-seq\n")
ggsave(file=makeOutFile("mRNA-cellExpr-dendroPolyA.pdf"),height=5,width=8)







#pc.pam.pam.melt1<-transform( pc.pam.pam.melt, new = ordered(.id, levels = names( sort(-table(value)))))
#pc.pam.pam.melt1 <- within(pc.pam.pam.melt, Position <- factor(.id, levels=names(sort(table(.id),decreasing=TRUE))))

ggplot(cd.melt.top1.df,aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~variable,ncol=2)+
  ggtitle("cshl transcript mRNA expression data\n(Top Expressed Transcript For Each Gene)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("mRNA-topTrans-cell-type-Expr.pdf"),height=24,width=16)

threeTimeTable <- as.data.frame(table(cd.top1.df$threeTimesRest))
ggplot(threeTimeTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
  theme_bw()+
  xlab("Frequency of Cell-type with 3x expression than rest")+
  ylab("cell-type and extraction method")+
  ggtitle("Number of transcripts in a cell-type with 3x RPKM values\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-threeTimeRest-countsPerCell.pdf"))

threeTimeTable <- as.data.frame(table(subset(cd.melt.top1.df, threeTimesRest == variable)[c("threeTimesRest")]))
ggplot(threeTimeTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
  theme_bw()+
  xlab("Frequency of Cell-type with 3x expression than rest")+
  ylab("cell-type and extraction method")+
  ggtitle(paste("Number of transcripts in a cell-type with 3x RPKM values\n(Top Expressed Transcript For Each Gene)\nN=",
                dim(subset(cd.melt.top1.df, threeTimesRest == variable)[c("threeTimesRest")])[1],
                "/",
                dim(cd.top1.df)[1],
                "(total)" ))
ggsave(file=makeOutFile("mRNA-topTrans-threeTimeRest-countsPerCell-noNone.pdf"))


ggplot(cd.melt.top1.df,aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  ggtitle("cshl transcript mRNA expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("mRNA-topTrans-threeTimeRest-exprInCellTypes.pdf"),height=24,width=16)

ggplot(cd.melt.top1.df[which(cd.melt.top1.df$threeTimesRest != "none"),],aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  ggtitle("cshl transcript mRNA expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("mRNA-topTrans-threeTimeRest-exprInCellTypes-NoNone.pdf"),height=24,width=16)

ggplot(cd.melt.top1.df[which(cd.melt.top1.df$threeTimesRest != "none"),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~threeTimesRest,ncol=2)+
  ggtitle("cshl mRNA transcript expression data\nfacetby three time rest expression\ncolor by threeTimesRest cell-type expr\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-facetThreeTimesRest.pdf"),height=24,width=16)

ggplot(cd.melt.top1.df[which(cd.melt.top1.df$tissSpec > 0.8),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("cshl mRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.8 tissue specifity\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-tissSpecGt0_8-facetMaxExpr.pdf"),height=24,width=16)


ggplot(cd.melt.top1.df[which(cd.melt.top1.df$tissSpec > 0.6),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("cshl mRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-tissSpecGt0_6-facetMaxExpr.pdf"),height=24,width=16)



#polar
ggplot(cd.melt.top1.df[which(cd.melt.top1.df$tissSpec > 0.6),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+coord_polar()+
  theme_bw()+
  ggtitle("cshl mRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-tissSpecGt0_6-combined-polar.pdf"),height=16,width=16)


#ggplot(cd.melt.top1.df,aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
#  geom_line()+coord_polar()+
 # theme_bw()+
#  ggtitle("cshl mRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity\n(Top Expressed Transcript For Each Gene)")
#ggsave(file=makeOutFile("mRNA-topTrans-combined-polar.pdf"),height=16,width=16)

ggplot(subset(cd.melt.top1.df,threeTimesRest != "none"),aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+coord_polar()+
  theme_bw()+
  ggtitle("cshl lncRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\nonly transcripts w/ > 0.6 tissue specifity\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-threeTimesRest-polar.pdf"),height=16,width=16)


#end polar



ggplot(cd.melt.top1.df,aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesRest)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("cshl mRNA transcript expression data\nfacet by sample with maximum expression\ncolor by threeTimesRest cell-type expr\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-facetMaxExpr.pdf"),height=24,width=16)

cd.melt.top1.df$cellTypeSpecific <- "no"
cd.melt.top1.df[which(cd.melt.top1.df$threeTimesRest == cd.melt.top1.df$variable),]$cellTypeSpecific <- "yes"

ggplot(cd.melt.top1.df[which(cd.melt.top1.df$value > 0),],aes(log(value),color=cellTypeSpecific))+geom_density()+
  theme_bw()+
  xlab("log RPMK")+
  ylab("density")+
  ggtitle("Cell type-specific expression\nComparing mRNA in cell type where expression is three time rest\nversus all other instances of mRNA expression\nwithin all cell type")
ggsave(file=makeOutFile("mRNA-cellTypeExpr-Comp-nonCellType.pdf"))

cd.melt.top1.df$ctsCat <- "none"
cd.melt.top1.df[which(cd.melt.top1.df$threeTimesRest == cd.melt.top1.df$variable),]$ctsCat <- "cell-type-this"
cd.melt.top1.df[intersect(which(cd.melt.top1.df$threeTimesRest != cd.melt.top1.df$variable), which(cd.melt.top1.df$threeTimesRest != "none")),]$ctsCat <- "cell-type-other"

ggplot(cd.melt.top1.df[which(cd.melt.top1.df$value>0),],aes(x=log(value),fill=ctsCat))+geom_density(alpha=I(0.4))+
  facet_wrap(~variable,ncol=2,scales="free")+
  ggtitle("Cell Type Expression:(top transcripts per gene)\nCell-type specificity of transcript is filled")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorBySpecificity.pdf"),height=24,width=12)


ggplot(cd.melt.top1.df[which(cd.melt.top1.df$value>0),],aes(x=log(value),colour=ctsCat))+geom_freqpoly(size=1,binwidth=0.5)+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("Cell Type Expression:(top transcripts per gene)\nCell-type specificity of transcript is colored")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorBySpecificity-freqPoly.pdf"),height=24,width=12)

#edit start

cd.melt.top1.df$maxExprCat <- "none"
cd.melt.top1.df[which(cd.melt.top1.df$maxExprType == cd.melt.top1.df$variable),]$maxExprCat <- "cell-type-this"
cd.melt.top1.df[intersect(which(cd.melt.top1.df$maxExprType != cd.melt.top1.df$variable), which(cd.melt.top1.df$maxExprType != "none")),]$maxExprCat <- "cell-type-other"

ctsCount.df <- ddply(subset(cd.melt.top1.df, value > 0),.(variable,maxExprCat),function(df)dim(df)[1])
colnames(ctsCount.df) <- c("RnaSeqExpr","maxExprCat","transcripts")
ctsCount.df$maxExprCat <- factor(ctsCount.df$maxExprCat,levels=c("cell-type-this","cell-type-other","none"))
ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=maxExprCat,colour=maxExprCat))+geom_histogram()+coord_flip()+ theme_bw()+
  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
  ggtitle("Cell-type Max Expressed Transcript\nFor mRNA transcripts found in {cell-type/polyA pulldown} pair")
ggsave(file=makeOutFile("mRNA-cellExpr-maxExprCountHistogram.pdf"),height=7,width=7)

cd.melt.top1.df$maxSpecCat <- "none"
cd.melt.top1.df[which(cd.melt.top1.df$threeTimesRest == cd.melt.top1.df$variable),]$maxSpecCat <- "cell-type-this"
cd.melt.top1.df[intersect(which(cd.melt.top1.df$threeTimesRest != cd.melt.top1.df$variable), which(cd.melt.top1.df$threeTimesRest != "none")),]$maxSpecCat <- "cell-type-other"

ctsCount.df <- ddply(subset(cd.melt.top1.df, value > 0),.(variable,maxSpecCat),function(df)dim(df)[1])
colnames(ctsCount.df) <- c("RnaSeqExpr","maxSpecCat","transcripts")
ctsCount.df$maxExprCat <- factor(ctsCount.df$maxSpecCat,levels=c("cell-type-this","cell-type-other","none"))
ggplot(ctsCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=maxSpecCat,colour=maxSpecCat))+geom_histogram()+coord_flip()+ theme_bw()+
  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
ggtitle("Cell-type Specific(ThreeTimesRest) Expressed Transcript\nFor mRNA transcripts found in {cell-type/polyA pulldown} pair")

ggsave(file=makeOutFile("mRNA-cellExpr-specifityCountHistogram.pdf"),height=7,width=7)
#edit end

##### JSD

cd.melt.top1.df$jsdCat <- "none"
cd.melt.top1.df[intersect(which(cd.melt.top1.df$maxExprType == cd.melt.top1.df$variable),which(cd.melt.top1.df$tissSpec > 0.6)),]$jsdCat <- "cell-type-this"
cd.melt.top1.df[intersect(which(cd.melt.top1.df$maxExprType != cd.melt.top1.df$variable),which(cd.melt.top1.df$tissSpec > 0.6)),]$jsdCat <- "cell-type-other"


ggplot(cd.melt.top1.df[which(cd.melt.top1.df$jsdCat != "none"),],aes(x=variable,y=log2(value),group=transcript_id,color=maxExprType)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExprType,ncol=2)+
  ggtitle("mNA\ncshl lncRNA top transcript expression data\nfacet and colory by JSD(0.6) Specificity metric\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-facetJSD.pdf"),height=24,width=16)





ggplot(cd.melt.top1.df[which(cd.melt.top1.df$value>0),],aes(x=log(value),fill=jsdCat))+geom_density(alpha=I(0.4))+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("mRNA\nCell Type Expression:(top transcripts per gene)\nfilled by category of transcript cell-type JSD specificity(0.6) of")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorByJsdSpecificity.pdf"),height=24,width=12)


ggplot(cd.melt.top1.df[which(cd.melt.top1.df$value>0),],aes(x=log(value),colour=jsdCat))+geom_freqpoly(size=1,binwidth=0.5)+
  facet_wrap(~variable,ncol=2,scale="free")+
  ggtitle("mRNA\nCell Type Expression:(top transcripts per gene)\nCell-type JSD specificity(0.6) of transcript is colored")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorByJSDSpecificity-freqPoly.pdf"),height=24,width=12)


jsdCount.df <- ddply(subset(cd.melt.top1.df, value > 0),.(variable,jsdCat),function(df)dim(df)[1])
colnames(jsdCount.df) <- c("RnaSeqExpr","jsdCat","transcripts")
jsdCount.df$jsdCat <- factor(jsdCount.df$jsdCat,levels=c("cell-type-this","cell-type-other","none"))
ggplot(jsdCount.df,aes(x=RnaSeqExpr,y=transcripts,fill=jsdCat,colour=jsdCat))+geom_histogram()+coord_flip()+ theme_bw()+
  scale_fill_manual(values = c("none" = brewer.pal(9,"Blues")[4],"cell-type-this" = brewer.pal(9,"Greys")[7],"cell-type-other" = brewer.pal(9,"Reds")[4]))+ 
  scale_colour_manual(values = c("none" = brewer.pal(9,"Blues")[5],"cell-type-this" = brewer.pal(9,"Greys")[8],"cell-type-other" = brewer.pal(9,"Reds")[5]))+
  ggtitle("mRNA\nCell-type JSD(0.6)\nFor transcripts found in {cell-type/polyA pulldown} pairs")
ggsave(file=makeOutFile("mRNA-cellExpr-jsdSpecCountHistogram.pdf"),height=7,width=7)

###END JS
####### JSD END



ggplot(cd.melt.top1.df[which(cd.melt.top1.df$value>0),],aes(x=log(value),fill=maxExprCat))+geom_density(alpha=I(0.4))+
  facet_wrap(~variable,ncol=2)+
  ggtitle("Cell Type Expression:(top transcripts per gene)\nCell-type max expression of transcript is filled")+
  theme_bw()
ggsave(file=makeOutFile("mRNA-cellExpr-colorByMaxExpr.pdf"),height=24,width=12)

maxExprTable <- as.data.frame(table(cd.top1.df$maxExprType))
ggplot(maxExprTable,aes(x=Var1,y=Freq))+geom_bar()+coord_flip()+
  theme_bw()+
  xlab("Frequency of Cell-type max expression of a lncRNA")+
  ylab("cell-type and extraction method")+
  ggtitle("Frequency of cell types transcript is max expressed in\n(Top Expressed Transcript For Each Gene)")
ggsave(file=makeOutFile("mRNA-topTrans-maxExpr-countsPerCell.pdf"),height=7,width=7)

#cdf 
cd.top1.df$isTissSpec <- "notThreeTimes"
cd.top1.df[which(cd.top1.df$threeTimesRest != "none"),]$isTissSpec <- "threeTimes"
tmp <- cd.top1.df
tmp$isTissSpec <- "combinedDistro"
tissSpec.df <- rbind(cd.top1.df,tmp)
#subset(lnc.top1.df, threeTimesRest != "none") -> lnc.spec.df
mydf_m <- ddply(tissSpec.df, .(isTissSpec), transform, ecd = ecdf(tissSpec)(tissSpec))
#mydf_m <- ddply(lnc.spec.df, .(variable), transform, ecd = ecdf(tissSpec)(tissSpec))

ggplot(mydf_m,aes(x = tissSpec, y = ecd)) + 
  geom_line(aes(group=isTissSpec,colour=isTissSpec,size=3)) +
  theme_bw()+
  ggtitle("CDF of tissSpecifity Measure(JS divergence)\nCummulative Frequency of ThreeTimesRest \nvs. Non-ThreeTimes rest detected(mRNA)")
ggsave(file=makeOutFile("mRNA-tissSpec-CDF.pdf"),height=7,width=7)


if(!("nTimesRest" %in% colnames(cd.top1.df))){
  cd.top1.df <- applyAndAppendFnDf(df=cd.top1.df,cols=expr.cols,FUN=function(x)nTimesRestSpecifity(x,expr.cols),"nTimesRest")
}

ggplot(cd.top1.df, aes(x=(tissSpec),y=log(nTimesRest)))+geom_point(alpha=I(0.4),size=1)+
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
  ggtitle(paste("Comparison of tissue specificity \nand (top expression n times greater than next)\nN=",dim(cd.top1.df)[1]))
ggsave(file=makeOutFile("mRNA-tissSpec-vs-nTimesRest.pdf"),height=7,width=7)


#############################################################################################################################################################################################
#
#    combined lncRNA/mRNA
#
#
##############################################################################################################################################################################################################

###########################
#     CELL/PULLDOWN
#     EXPR 
########################### ggplot(lnc.comb[which(lnc.comb$expression > 0),],aes(x=log(expression),color=seqPullDown))+geom_freqpoly(size=1,binwidth=0.5)+

lnc.melt.top1.df$transType = "lncRNA"
cd.melt.top1.df$transType = "mRNA"
if (all.equal(colnames(lnc.melt.top1.df),colnames(cd.melt.top1.df))){combined.melt.top1.df <- rbind(lnc.melt.top1.df,cd.melt.top1.df)}

ggplot(combined.melt.top1.df,aes(x=log2(value),color=transType))+
  geom_freqpoly(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~variable,ncol=2)+
  ggtitle("cshl transcript lncRNA vs. mRNA expression data\n(Top Expressed Transcript For Each Gene)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("combined-topTrans-cell-type-Expr.pdf"),height=24,width=12)

ggplot(combined.melt.top1.df,aes(x=log2(value),color=transType))+
  geom_freqpoly(binwidth=0.4)+
  theme_bw()+
  ggtitle("cshl transcript lncRNA vs. mRNA expression data\n(Top Expressed Transcript For Each Gene)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=makeOutFile("combined-topTrans-Expr.pdf"),height=7,width=7)

as.data.frame(as.list(diag(cor(lnc.lnpa.df[expr.cells],lnc.lpa.df[expr.cells])))) -> lnc.cor
lnc.cor$type = "lncRNA"
as.data.frame(as.list(diag(cor(cd.lnpa.df[expr.cells],cd.lpa.df[expr.cells])))) -> cd.cor
cd.cor$type = "mRNA"
comb.cor <- melt(rbind(cd.cor, lnc.cor),id.var="type")
colnames(comb.cor) <- c("type","cell","correlationCoeff")
ggplot(comb.cor,aes(x=cell, y=correlationCoeff))+geom_histogram()+facet_wrap(~type,nrow=2)+theme_bw()+coord_flip()+ggtitle("Correlation Coeff between polyA and nonPolyA\nfor each cell type,transcript type")
ggsave(file=makeOutFile("combined-correlation-exprVals.pdf"),height=7,width=7)


cd.ecdf <- ecdf(cd.melt.top1.df$value)
with(cd.top1.df,ecdf(get(expr.cols[1]))) -> col1
qq.col1 <- data.frame(x=as.vector(sapply(cd.top1.df[expr.cols[1]],col1)),y=as.vector(sapply(cd.top1.df[expr.cols[1]],cd.ecdf)))
ggplot(qq.col1,aes(x=x,y=y))+geom_point()+geom_abline(slope=1,intercept=0)+xlim(0,1)+ylim(0,1)
ggsave(makeTmpOutFile("qqplotTest.pdf"),height=7,width=7)

