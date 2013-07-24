home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()
library(hexbin)
library(ggplot2)
library(plyr)
library(reshape)
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
outdir <- getFullPath("plots/cshlRnaSeq")

derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"))
lnc.gene.trans <- derrien.df[c("LncRNA_GeneId","LncRNA_Txid")]
cshl.trans.df <- read.csv(file=getFullPath("data/cshlLncTransWithStats.tab"),sep=" ")
cshl.trans.df <- cshl.trans.df[which(cshl.trans.df$average > 0),]
cshl.trans.df <- merge(cshl.trans.df,lnc.gene.trans,by.x="transcript_id",by.y="LncRNA_Txid")


csh1.gene.var.df <- ddply(cshl.trans.df,.(LncRNA_GeneId),numcolwise(var))
csh1.gene.mean.df <- ddply(cshl.trans.df,.(LncRNA_GeneId),numcolwise(mean))

csh1.gene.var.melt.df <- melt(csh1.gene.var.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562","LncRNA_GeneId")],id="LncRNA_GeneId")
csh1.gene.mean.melt.df <- melt(csh1.gene.mean.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562","LncRNA_GeneId")],id="LncRNA_GeneId")

colnames(csh1.gene.mean.melt.df)<- c("LncRNA_GeneId","variable","meanValue")
csh1.gene.var.melt.df$meanValue = csh1.gene.mean.melt.df$meanValue # add the mean to the variance data.frame

genes<-cshl.trans.df$LncRNA_GeneId
cshl.trans.df$dummyGene <- sample(genes,length(genes),replace=FALSE)
csh1.gene.var.dummy.df <- ddply(cshl.trans.df,.(dummyGene),numcolwise(var))
csh1.gene.mean.dummy.df <- ddply(cshl.trans.df,.(dummyGene),numcolwise(mean))
csh1.gene.var.melt.dummy.df <- melt(csh1.gene.var.dummy.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562","dummyGene")],id="dummyGene")
csh1.gene.mean.melt.dummy.df <- melt(csh1.gene.mean.dummy.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562","dummyGene")],id="dummyGene")
colnames(csh1.gene.mean.melt.dummy.df)<- c("dummyGene","variable","meanValue")
csh1.gene.var.melt.dummy.df$meanValue = csh1.gene.mean.melt.dummy.df$meanValue

cshl.trans.df$variance <- apply(cshl.trans.df[c("GM12878","HELAS3", "K562","HEPG2","H1HESC")],1,function(...){var(c(...))})
cshl.trans.df$stddev <- apply(cshl.trans.df[c("GM12878","HELAS3", "K562","HEPG2","H1HESC")],1,function(...){sqrt(var(c(...)))})

mdata <- melt(cshl.trans.df[c("transcript_id","GM12878","HELAS3", "K562","HEPG2","H1HESC","entropy","threeTimesSpecifity","maxExpr","tissSpec")], id=c("transcript_id","entropy","threeTimesSpecifity","maxExpr","tissSpec"))
data$quartileE <- with(mdata, cut(entropy, breaks=quantile(entropy, probs=seq(0,1, by=0.33)), include.lowest=TRUE))
mdata$quartileE <- with(mdata, cut(entropy, breaks=quantile(entropy, probs=seq(0,1, by=0.33)), include.lowest=TRUE))

mdata$quartileTS <- with(mdata, cut(tissSpec, breaks=quantile(tissSpec, probs=seq(0,1, by=0.33)), include.lowest=TRUE))


nozero <- cshl.trans.df[which( cshl.trans.df$max == 0),]
cshl.trans.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562","LncRNA_GeneId")]
split(cshl.trans.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562")],cshl.trans.df$LncRNA_GeneId)
cor(cshl.trans.df[which(cshl.trans.df$LncRNA_GeneId == "ENSG00000256995.1"),c("HELAS3","GM12878","H1HESC","HEPG2","K562")])
#ENSG00000229404.1
#ENSG00000229401.1
aveDistance <- function(df){
  dist(scale(df))
}

dff <- function(df){cm <- colMeans(df[-1])}
ef <- function(df)ifelse(dim(df)[1]>2,eigen(cor(df))$val[1],-1)
d <- function(df)ave(dist(df))[1]
linD <- function(df){apply(df[-1],1,function(x){ave(sqrt(sum(x * colMeans(df[-1]))))})}
linD1 <- function(df){as.data.frame(apply(df[-1],1,function(x){ave(sqrt(sum(x * colMeans(df[-1]))))}))}
sumDist <- function(df)as.list(dist(df))

ldply(split(cshl.trans.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562")],cshl.trans.df$LncRNA_GeneId),linD1)->x
colnames(x)<-c("LncRNA_GeneId","distance")
as.data.frame(as.table(x))->check
y<-merge(x,cshl.trans.df)
ggplot(y,aes(log(distance)))+geom_bar()


sapply(split(cshl.trans.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562")],cshl.trans.df$LncRNA_GeneId),sumDist)->x.forDistance
ldply(x.forDistance,function(li){data.frame(distance=unlist(li))})->x.dist.df
ggplot(x.dist.df,aes(x=log(distance)))+geom_bar()

#cshl.trans.melt.df <- merge(mdata,cshl.trans.df[c("transcript_id")])

# stdev VS. mean
ggplot(csh1.gene.var.melt.df[which(!is.na(csh1.gene.var.melt.df$value)),],aes(x=log(meanValue),y=log(sqrt(value))))+geom_point()+facet_wrap(~variable)+
  ggtitle("cshl transcript isoform lncRNA expression mean vs. standard deviation\n(calculated over all the transcripts of one gene")+
  ylab("log2(standard deviation of each transcript RPKM for each gene,whole cell)")+
  xlab("log2(mean of each transcript RPKM for each gene,whole cell)")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Gene-Stdev-Mean.pdf"))

genes<-csh1.gene.var.melt.df$LncRNA_GeneId
csh1.gene.var.melt.df$dummyGene <- sample(genes,length(genes),replace=FALSE)

ggplot(csh1.gene.var.melt.dummy.df[which(!is.na(csh1.gene.var.melt.dummy.df$value)),],aes(x=log(meanValue),y=log(sqrt(value))))+geom_point()+facet_wrap(~variable)+
  ggtitle("DUMMY: cshl transcript isoform lncRNA expression mean vs. standard deviation\n(Genes have be shuffled")+
  ylab("log2(standard deviation of each transcript RPKM for each gene,whole cell)")+
  xlab("log2(mean of each transcript RPKM for each gene,whole cell)")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Gene-Stdev-Mean-DUMMY.pdf"))

sim.df <- data.frame(mean=rep(0,3000),stdev=rep(0,3000))
colnames(sim.df) <-c("mean","stdev")
for(i in 1:3000){
  x <- 2^rexp(10,1/4)
  sim.df$stdev[i] <- sqrt(var(x))
  sim.df$mean[i] <- mean(x)
}
ggplot(sim.df,aes(x=log2(mean),y=log2(stdev)))+geom_point()
ggsave("/Users/adam/Desktop/sim.pdf")

ggplot(csh1.gene.var.melt.df,aes(x=log(value)))+geom_bar()+facet_wrap(~variable)+
ggtitle("cshl transcript lncRNA expression variance\n")+
  xlab("log2(variance of each transcript RPKM for each gene,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Gene-Variance.pdf"))

ggplot(mdata,aes(x=log(value)))+geom_bar(binwidth=0.4)+facet_wrap(~variable)+
  ggtitle("cshl transcript lncRNA expression ")+
  xlab("log2(transcript RPKM for each gene,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Trans-Expr.pdf"))



ggplot(mdata,aes(x=log2(value)))+
  geom_bar(binwidth=0.4)+
  theme_bw()+
  facet_wrap(~threeTimesSpecifity)+
  ggtitle("cshl transcript lncRNA expression data\ngrouped by cell with transcript expression\n 3x greater than other cell-types")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Trans-cts.pdf"))

ggplot(mdata,aes(x=log2(value),fill=factor(threeTimesSpecifity )))+
  geom_bar(binwidth=0.4,position="fill")+
  theme_bw()+facet_wrap(~variable)+
  ggtitle("cshl lncRNA transcript expression data\nColor by transcripts' 3x expression cell-type")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Trans-cst2.pdf"))

ggplot(mdata,aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesSpecifity)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~quartileTS)+
  ggtitle("cshl lncRNA transcript expression data\nfacetby tissue specifity\ncolor by three times rest cell-type expr")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Trans-cst-TSmutualI.pdf"))

ggplot(mdata[which(mdata$tissSpec > 0.8),],aes(x=variable,y=log2(value),group=transcript_id,color=threeTimesSpecifity)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExpr)+
  ggtitle("cshl lncRNA transcript expression data\nfacetby tissue specifity\ncolor by three times rest cell-type expr\nonly transcripts w/ > 0.8 tissue specifity")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Trans-cst-TSmutual-high-specificity.pdf"))


mdata <-  mdata[which(mdata$maxExpr != "None"),]
ggplot(mdata[which(mdata$threeTimesSpecifity != "None"),],aes(x=variable,y=log2(value),group=transcript_id,color=tissSpec)) + 
  geom_line()+
  theme_bw()+
  facet_wrap(~maxExpr)+
  ggtitle("cshl lncRNA transcript expression data\nfacetby tissue specifity\ncolor by JS relative entropy in max cell line")+
  scale_colour_gradient(low="red", high="green")
ggsave(file=getFullPath("plots/cshl/cslh-lncRNA-Trans-threeTimesRestExpressionProfile.pdf"))

ggplot(cshl.trans.df,aes(x=log2(average),y=log2(variance)))+geom_point()
ggsave(file=getFullPath("plots/cshlRnaSeq/lncRNA-Expr-ave-vs-variance.pf"))

ggplot(as.data.frame(table(derrien.df$LncRNA_GeneId)),aes(x=Freq))+geom_histogram(binwidth=1)+
  ggtitle("Isoforms per lncRNA gene")+
  ylab("count")+
  xlab("isoforms per gene")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Isoforms-Per-Gene.pdf"))

as.data.frame(table(derrien.df$LncRNA_GeneId))->geneIsoform.df
colnames(geneIsoform.df)<-c("LncRNA_GeneId","isoforms")
geneTransIsoCounts.df <- merge(geneIsoform.df,derrien.df[c("LncRNA_GeneId", "LncRNA_Txid")])
mdata3 <- merge(mdata,geneTransIsoCounts.df,by.x="transcript_id",by.y="LncRNA_Txid")

mdata3$isoforms = factor(mdata3$isoforms)
mdata3[c("value","isoforms","variable")]->mdata4
ggplot(mdata4,aes(x=log2(value),fill=factor(isoforms)))+
    geom_bar(binwidth=1)+
    facet_wrap(~variable )+
    theme_bw()+
    ggtitle("cshl lncRNA transcript expression data\ncolor=number of isoforms")+
    xlab("log2(RPKM,whole cell)")+
    ylab("count")
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-cell-types.pdf"))

mdata5 <- mdata4[which(mdata4$isoforms != 1),]
ggplot(mdata5,aes(x=log2(value),fill=factor(isoforms)))+
  geom_bar(binwidth=1)+ ggtitle("cshl lncRNA transcript expression data\ncolor=number of isoforms\n(genes w/ > 1 isoform)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-cell-types-non-one.pdf"))

ggplot(mdata5,aes(x=log2(value),fill=factor(isoforms)))+
  geom_bar(binwidth=1,position="fill")+ ggtitle("cshl lncRNA gene expression data\ncolor=number of isoforms\n(genes w/ > 1 isoform)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-cell-types-non-one-fill.pdf"))


ggplot(mdata5,aes(x=value,fill=factor(isoforms)))+
  geom_bar(binwidth=0.1,position="fill")+ ggtitle("cshl lncRNA transcription expression data\ncolor=number of isoforms\n(genes w/ > 1 isoform)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")+
  xlim(1,4)
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-cell-types-non-one-range.pdf"))

ggplot(mdata4,aes(x=log2(value),fill=factor(isoforms)))+
  geom_bar(binwidth=1)+ ggtitle("cshl lncRNA transcript expression data\ncolor=number of isoforms\n(genes w/ >= 1 isoform)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-Allcell-types-.pdf"))

ggplot(mdata4,aes(y=log2(value),x=factor(isoforms)))+
  geom_boxplot()+ ggtitle("cshl lncRNA transcript expression data\ncolor=number of isoforms\n(genes w/ >= 1 isoform)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")+
  theme_bw()
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-Allcell-types-boxplot.pdf"))


ggplot(mdata4,aes(x=log2(value),fill=factor(isoforms)))+
  geom_bar(binwidth=1,position="fill")+ ggtitle("cshl lncRNA transcript expression data\ncolor=number of isoforms\n(genes w/ >= 1 isoform)")+
  xlab("log2(RPKM,whole cell)")+
  ylab("count")
ggsave(file=getFullPath("plots/derrien2012/lnc-RNA-isoformPerGene-expr-in-Allcell-types-fill.pdf"))
  

count<-sapply(split(geneIsoform.df,geneIsoform.df$isoforms),function(df)sum(df$isoforms))
data.frame(count=as.vector(count),isoforms=seq(1,length(as.vector(count)),1 ))-> transcriptGene.df
ddply(transcriptGene.df,.(isoforms),function(df)data.frame(x=rep(df$isoform,df$count),y=rep(1,df$count)))->transPerGeneByTrans.df

ggplot(transPerGeneByTrans.df,aes(x=x)) +
  geom_bar(binwidth=1)+ theme_bw()+
  ggtitle("lncRNA transcript's gene isoform count")+
  xlab("transcript's gene isoform count")+
  ylab("count")
ggsave(file=getFullPath("plots/cshl/lnc-RNA-isoformPerGene-expr-combined.pdf"))

ggplot(transPerGeneByTrans.df,aes(x=x))+
  stat_ecdf()+
  theme_bw()+
  xlab("isoforms within gene transcript originates from")+
  ylab("cummulative probability")+
  ggtitle("lncRNA transcript's gene number of isoform")
ggsave(file=getFullPath("plots/cshl/lnc-RNA-isoformPerGene-cdf.pdf"))



princomp(cshl.trans.df[c("HELAS3","GM12878","H1HESC","HEPG2","K562")])->lnc.pca
as.data.frame(lnc.pca$scores) -> lnc.pca.factors
lnc.pca.factors$transcript_id <- cshl.trans.df$transcript_id
lnc.pca.melt <- melt(lnc.pca.factors)
ggplot(lnc.pca.melt,aes(x=value,fill=variable))+geom_density()
ggsave("~/Desktop/pcacolors.pdf")

ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2)))+stat_binhex()
ggplot(lnc.pca.melt,aes(x=value,fill=variable))+theme_bw()+geom_density(alpha=I(0.4))+xlim(-0.1,0.1)
ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2)))+geom_point()+
theme_bw()+
  xlab("Component 1")+
  ylab("Component 2")+
  ggtitle("lncRNA transcript expression PCA project on two components")
ggsave(file=getFullPath("plots/cshl/lnc-RNA-expression-PCA.pdf"))

