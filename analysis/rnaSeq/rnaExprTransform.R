home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }


file = getFullPath("data/combExp.tab")
file1 = getFullPath("data/combinedExprLncRNA.txt")
lnc.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"))
df <- read.table(file=file1,  header=TRUE, stringsAsFactors=FALSE)
df.1 <- df
cells = c("HELAS3","GM12878","H1HESC","HEPG2","K562")

entropy <- function(...){
  x = c(...)
  x = x / sum(x)
  y = log2(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}
entropyVec <- function(x){
  x = x / sum(x)
  y = log2(x) * x
  y[is.na(y)]=0
  sum(y) * -1
}
JS = function(x,y){
  x = x/ sum(x)
  y = y / sum(y)
  a = (x + y)/2
  entropyVec(a )- ((entropyVec(x)-entropyVec(y))/2) -> z
  z}

JSsp <- function(e1,e2){
  1 - sqrt(JS(e1,e2))
}

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




df.1$entropy = apply(df[,-1],1,entropy)
df.1$variance = apply(cshl.trans.df[c("GM12878","HELAS3", "K562","HEPG2","H1HESC")],1,function(...){var(c(...))})
df.1$min = apply(df[,-1],1,min)
df.1$max = apply(df[,-1],1,max)
fmean<-function(...){mean(c(...))}
df.1$average = apply(df[,-1],1,fmean)
maxExpr <- function(...){
  ifelse( sum(c(...)) > 0.00000000001,
  c("HELAS3","GM12878","H1HESC","HEPG2","K562")[which(c(...) == max(c(...)))],
  "None"
  )
}
df.1$maxExpr = apply(df[,-1],1,maxExpr)
#df.2<-df.1[which(df.1$maxExpr != "None"),]
df.2 <- df.1

notExpressed <- "notExpressed" # change this to NA? or FALSE???
df.2$HELAS3.Expr = as.character(cut(df.2$HELAS3, quantile(df.2[which(df.2$HELAS3 > 0),]$HELAS3,(0:3)/3),labels=c("low","mid","high")))
df.2[which(is.na(df.2$HELAS3.Expr)),]$HELAS3.Expr <- notExpressed

df.2$GM12878.Expr = as.character(cut(df.2$GM12878, quantile(df.2[which(df.2$GM12878 > 0),]$GM12878,(0:3)/3),labels=c("low","mid","high")))
df.2[which(is.na(df.2$GM12878.Expr)),]$GM12878.Expr <- notExpressed

df.2$H1HESC.Expr = as.character(cut(df.2$H1HESC, quantile(df.2[which(df.2$H1HESC > 0),]$H1HESC,(0:3)/3),labels=c("low","mid","high")))
df.2[which(is.na(df.2$H1HESC.Expr)),]$H1HESC.Expr <- notExpressed

df.2$HEPG2.Expr = as.character(cut(df.2$HEPG2, quantile(df.2[which(df.2$HEPG2 > 0),]$HEPG2,(0:3)/3),labels=c("low","mid","high")))
df.2[which(is.na(df.2$HEPG2.Expr)),]$HEPG2.Expr <- notExpressed

df.2$K562.Expr = as.character(cut(df.2$K562, quantile(df.2[which(df.2$K562 > 0),]$K562,(0:3)/3),labels=c("low","mid","high")))
df.2[which(is.na(df.2$K562.Expr)),]$K562.Expr <- notExpressed

threeTimesRestSpecifity <- function(...){
  x <- c(...)
  ifelse( sort(x,decreasing=TRUE)[1] > (3* sort(x,decreasing=TRUE)[2]),
      c("HELAS3","GM12878","H1HESC","HEPG2","K562")[which(x == max(x))],
          "None")}
          
df.2$threeTimesSpecifity = apply(df.2[2:6],1,threeTimesRestSpecifity)
df.2$tissSpec = apply(df.2[2:6],1,tissSpec)

df.2$expressed = "no"
df.2[which(df.2$maxExpr != "None"),]$expressed <- "yes"
merged.df<-merge(lnc.df,df.2,by.x="LncRNA_Txid",by.y="transcript_id")

exportAsTable(df.2,getFullPath("/Users/adam/Desktop/combinedLncExprWithStats.tab"))
exportAsTable(merged.df,getFullPath("/Users/adam/Desktop/Gencode_lncRNAsv7_summaryTable_05_02_2012_lncExression.tab"))

geneIsoform.df$isoformSqr<-apply(geneIsoform.df[2],1,function(x){x*x})
geneI.df <- ddply(geneIsoform.df,.(isoforms),subset,sum)
ggplot(geneIsoform.df,aes(x=isoformSqr))+geom_bar()
