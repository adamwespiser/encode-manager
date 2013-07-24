home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }

df <- read.csv(file=getFullPath("data/gencodeV7lncRNAexpr.csv"))
lnc.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"))
 
library(ggplot2)
library(plyr)
library(reshape)

information <- function(x){
  x = x / sum(x)
  y = log2(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}



JS = function(x,y){
  x = x/ sum(x)
  y = y / sum(y)
  a = (x + y)/2
  information(a )- ((information(x)-information(y))/2) -> z
  z}


mdata <- melt(df, id=c("filename"))
ggplot(mdata,aes(x=value))+geom_bar(binwidth=0.1)+xlim(0.1,4)+theme_bw()+facet_wrap(~variable)+ggtitle("derrien2012 lncRNA gene expression data")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr.pdf"))

df.2 <- df[c("filename", "HEPG2", "GM12878", "K562", "HELAS3", "H1HESC")]

#df.2 <-transform(df.2, i =information(c(NHEK,GM12878,K562,HELAS3,H1HESC)))

df.2$sum = df.2$HEPG2 + df.2$GM12878 + df.2$K562  +  df.2$HELAS3 +  df.2$H1HESC

df.2$HEPG2p = df.2$HEPG2/df.2$sum 
df.2$GM12878p = df.2$GM12878/df.2$sum  
  df.2$K562p =   df.2$K562/df.2$sum  
  df.2$HELAS3p = df.2$HELAS3/df.2$sum  
  df.2$H1HESCp =  df.2$H1HESC/df.2$sum
df.2 <- df.2[which(df.2$sum > 0),]

df.2 = transform(df.2, iscore = -1 * (  ifelse(is.na(log2(df.2$HEPG2p)*df.2$HEPG2p), 0, log2(df.2$HEPG2p)*df.2$HEPG2p) +   
                        ifelse(is.na(log2(df.2$GM12878p)*df.2$GM12878p), 0, log2(df.2$GM12878p)*df.2$GM12878p) +
                        ifelse(is.na(log2(df.2$K562p)*df.2$K562p), 0, log2(df.2$K562p)*df.2$K562p)+
                        ifelse(is.na(log2(df.2$HELAS3p)*df.2$HELAS3p), 0, log2(df.2$HELAS3p)*df.2$HELAS3p)+
                        ifelse(is.na(log2(df.2$H1HESCp)*df.2$H1HESCp),0,log2(df.2$H1HESCp)*df.2$H1HESCp)))
# df.2 = transform(df.2, i = )
 

df.3 <- df.2[c("filename","GM12878","HELAS3", "K562","HEPG2","H1HESC", "iscore")]
mdata2 <- melt(df.3, id=c("filename","iscore"))
mdata2$quartile <- with(mdata2, cut(iscore, breaks=quantile(iscore, probs=seq(0,1, by=0.25)), include.lowest=TRUE))

 ggplot(mdata2,aes(x=log2(value),fill=quartile))+
   geom_bar(binwidth=0.4)+
   theme_bw()+
   facet_wrap(~variable)+
   ggtitle("derrien2012 gene lncRNA expression data\nquartile=shanon entropy")+
   xlab("log2(RPKM,whole cell)")+
   ylab("count")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr-SE.pdf"))

 ggplot(mdata2,aes(x=log2(value),fill=quartile))+
   geom_bar(binwidth=0.4,position="fill")+
   theme_bw()+facet_wrap(~variable)+
   ggtitle("derrien2012 lncRNA gene expression data\nquartile=shanon entropy")+
   xlab("log2(RPKM,whole cell)")+
   ylab("count")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr-SEstack.pdf"))

 ggplot(mdata2,aes(x=log2(value))) + 
  geom_bar(binwidth=0.2)+
  theme_bw()+
   facet_wrap(~variable)+
   ggtitle("derrien2012 lncRNA gene expression data")+ xlab("log2(RNA expr)")+
   xlab("log2(RPKM,whole cell)")+
   ylab("count")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr2.pdf"))


 ggplot(mdata2,aes(x=variable,y=value,color=quartile,size=5))+ 
   geom_jitter()+
   theme_bw() +
   ggtitle("derrien2012 lncRNA gene expression data")+ xlab("log2(RNA expr)")+
   ylab("RPKM,whole cell")+
   xlab("cell-type")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr-jitter.pdf")) 

 ggplot(mdata2,aes(x=variable,y=log2(value),group=filename,color=value)) + 
    geom_line()+
    theme_bw()+
    facet_wrap(~quartile)
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr-Cells.pdf")) 

 l.df <- merge(mdata2,lnc.df, by.x="filename", by.y="LncRNA_Txid")
 
 ggplot(l.df,aes(x=variable,y=log2(value),group=filename,color=phastcons_mammal_transcript_score))+
    geom_line()+
    theme_bw()+
    facet_wrap(~quartile)
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr-Cells-trans.pdf")) 

ggplot(l.df[which(l.df$phastcons_mammal_transcript_score > 300), ],aes(x=variable,y=log2(value),group=filename,color=phastcons_mammal_transcript_score))+
  geom_line()+theme_bw()+
  facet_wrap(~quartile)+
  scale_colour_gradient(low="red", high="green")
ggsave(file=getFullPath("plots/derrien2012/derrien2012-lncRNA-Gene-Expr-Cells-trans.pdf"))  
 
 
 num<-function(x){if(x=="NHEK"){ y<-4};if(x=="GM12878"){ y<-2};if(x=="K562"){ y<-3};if(x=="HELAS3"){y<-1};if(x=="H1HESC"){ y<-5};y}
 l.df$number <- 0
 l.df[which(l.df$variable == "NHEK"), ]$number <- 4
 l.df[which(l.df$variable == "GM12878"), ]$number <- 2
 l.df[which(l.df$variable == "K562"), ]$number <- 3
 l.df[which(l.df$variable == "HELAS3"), ]$number <- 1
 l.df[which(l.df$variable == "H1HESC"), ]$number <- 5
 
 expand.grid(x=seq(1,30),y=seq(1,30))->         dd
 transform(dd , sum= x + y)-> dd
 transform(dd, xp=x/sum, yp=y/sum)           -> dd
 transform(dd, xlog=log2(xp), ylog=log2(yp)) -> dd
 transform(dd, xv=xlog*xp, yv=ylog*yp)       -> dd
 transform(dd, info= -1*(xv+yv))             -> dd
 ggplot(dd,aes(x=x,y=y,size=info*.5,alpha=0.4,color=info))+geom_point()+theme_bw()
 
 