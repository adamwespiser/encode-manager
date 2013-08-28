#home <- Sys.getenv("HOME")
home <- "/Users/adam"

projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

#boiler plate helpers
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)



runNHEKCompare = function(outdir = getFullPath("plots/rnaSeq-QC/"), file = getFullPath("data/lncRnaExpr_ml.tab")){

df = readInTable(file)

cols = c("NHEK.longPolyA.1","NHEK.longPolyA","NHEK.longNonPolyA.1","NHEK.longNonPolyA")


#long Poly A 
df$NHEK.longPolyA.norm.1 = ifelse( df$NHEK.longPolyA.1 + 1 == 0, 0, log(df$NHEK.longPolyA.1 + 1))
df$NHEK.longPolyA.norm   = ifelse( df$NHEK.longPolyA   + 1 == 0, 0, log(df$NHEK.longPolyA   + 1))


model.lpa = lm(formula = df$NHEK.longPolyA.norm ~ df$NHEK.longPolyA.norm.1)

ggplot(df, aes(x=log(log(NHEK.longPolyA.1 + 1)),y=log(log(NHEK.longPolyA+ 1)),color=factor(label),size=(3 + label)))+geom_point() + theme_bw()+
  #ggplot(df, aes(x=log(NHEK.longPolyA.norm.1),y=log(NHEK.longPolyA.norm),color=factor(label),size=(3 + label)))+geom_point() + theme_bw()+
  
  geom_abline(slope= model.lpa$coefficients[2],intercept= model.lpa$coefficients[1]) +
  ggtitle(paste("RNA-seq biological replicates\nlong Poly A w/ psuedocount:(expr = expr + 1)\nCorrelation = ",round(cor( df$NHEK.longPolyA.1 , df$NHEK.longPolyA),digits=5)))
ggsave(file=paste(outdir,"NHEK-longPolyA-compare.pdf"))



#df = df[which(df$NHEK.longNonPolyA.norm.1 * df$NHEK.longNonPolyA.norm > 0), ]

df$NHEK.longNonPolyA.norm.1 =ifelse( df$NHEK.longNonPolyA.1 + 1 == 0, 0, log(df$NHEK.longNonPolyA.1 + 1))
df$NHEK.longNonPolyA.norm   =ifelse( df$NHEK.longNonPolyA   + 1 == 0, 0, log(df$NHEK.longNonPolyA   + 1))

model.lnpa = lm(df$NHEK.longNonPolyA.norm ~ df$NHEK.longNonPolyA.norm.1 )

ggplot(df, aes(x=log(log(NHEK.longNonPolyA.1+1)),y=log(log(NHEK.longNonPolyA+1)),color=factor(label),size=(3 + label)))+geom_point() + theme_bw()+
#ggplot(df, aes(x=log(NHEK.longNonPolyA.norm.1),y=log(NHEK.longNonPolyA.norm),color=factor(label),size=(3 + label*1))) + geom_point() + theme_bw()+
  geom_abline(slope= model.lnpa$coefficients[2],yintercept= model.lnpa$coefficients[1])+
 ggtitle(paste("RNA-seq biological replicates\nlong non-Poly A\nCorrelation = ",round(cor( df$NHEK.longNonPolyA.1 , df$NHEK.longNonPolyA),digits=5)))
ggsave(file=paste(outdir,"NHEK-longNonPolyA-compare.pdf"))

}


#df.cutoff = df[which(df$NHEK.longPolyA.1 > exp(-4) & df$NHEK.longPolyA > exp(-4)),]
#df.cutoff = df[which(df$NHEK.longPolyA.1 > 0 & df$NHEK.longPolyA > 0),]


plotTwoColsCompare = function(df, col1, col2, title, outdir){
d.tmp = data.frame(
c1 = df[[col1]],
c2 = df[[col2]],
c1log = log(df[[col1]] + 1),
c2log = log(df[[col2]] + 1),
c1loglog = log(log(df[[col1]] + 1) + 1),
c2loglog = log(log(df[[col2]] + 1) + 1),
label = df[["label"]])

col1f = gsub("\\.","",col1)
col2f = gsub("\\.","",col2)

colsCorr = round(cor(d.tmp$c1,d.tmp$c2),4)

ggplot(d.tmp, aes(x=c1,y=c2,color=factor(label), size=label + 1)) + geom_point() + theme_bw()+
	ggtitle(paste(title,"\ncorrelation: ",colsCorr,sep=""))+ 
  xlab(col1) + ylab(col2) 
ggsave(file=paste(outdir,"-compare.pdf",sep=""))


ggplot(d.tmp, aes(x=c1log,y=c2log,color=factor(label), size=label + 1)) + geom_point() + theme_bw()+
	ggtitle(paste(title,"\ncorrelation: ",colsCorr,sep=""))+ 
  xlab(paste("log(",col1,"+ 1" )) + ylab(paste("log(",col2,"+ 1)"))
ggsave(file=paste(outdir,"-log-compare.pdf",sep=""))

ggplot(d.tmp, aes(x=c1loglog,y=c2loglog,color=factor(label), size=label + 1)) + geom_point() + theme_bw()+
	ggtitle(paste(title,"\ncorrelation: ",colsCorr,sep=""))+ 
  xlab(paste("log(log(",col1,"+ 1) + 1)" )) + ylab(paste("log(log(",col2,"+ 1) + 1)"))
ggsave(file=paste(outdir,"-loglog-compare.pdf",sep=""))

d.tmp$c1only = ifelse(d.tmp$c1 > 0 & d.tmp$c2 == 0,1,0)
d.tmp$c2only = ifelse(d.tmp$c2 > 0 & d.tmp$c1 == 0,1,0)

print("heere")
ggplot(d.tmp[which(d.tmp$c1only == 1),], aes(c1)) + geom_density() + theme_bw()+
  ggtitle(paste(title,"\n",col1,"> 0,",col2,"=",0)) + xlab(col2)
ggsave(file=paste(outdir,col1f,"-only-density.pdf",sep=""))

ggplot(d.tmp[which(d.tmp$c2only == 1),], aes(c2)) + geom_density() + theme_bw()+
  ggtitle(paste(title,"\n",col2,"> 0,",col1,"=",0)) + xlab(col2)
ggsave(file=paste(outdir,col2f,"-only-density.pdf",sep=""))

}


funCompare2 = function(){
file.expr = getFullPath("data/lncRnaExpr_ml.tab")
df = readInTable(file.expr)

cols = c("NHEK.longPolyA.1","NHEK.longPolyA","NHEK.longNonPolyA.1","NHEK.longNonPolyA")

plotTwoColsCompare(df, col1=cols[1], col2=cols[2], title="NHEK longPolyA compare", outdir= getFullPath("plots/rnaSeq-QC/NHEK-longPolyA"))

plotTwoColsCompare(df, col1=cols[3], col2=cols[4], title="NHEK longNonPolyA compare", outdir= getFullPath("plots/rnaSeq-QC/NHEK-longNonPolyA"))

}



test = function(){
 d = data.frame(x = 1:1000/1000 + 100,y=1:1000+runif(1000)*100+10000)
  m = lm(formula=d$y~ d$x )
  ggplot(d,aes(x=x,y=y))+geom_point()+geom_abline(slope=m$coefficients[2],intercept= m$coefficients[1],size=4,color="green")+theme_bw()
  
}

