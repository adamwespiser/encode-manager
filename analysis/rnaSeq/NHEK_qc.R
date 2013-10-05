source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/createMLdataset.R")



runNHEKCompare = function(outdir = getFullPath("plots/rnaSeq-QC/"), file = getFullPath("data/lncRnaExpr_ml.tab")){

df = readInTable(file)

cols = c("NHEK.longPolyA.1","NHEK.longPolyA","NHEK.longNonPolyA.1","NHEK.longNonPolyA")


#long Poly A 
df$NHEK.longPolyA.norm.1 = ifelse( df$NHEK.longPolyA.1 + 1 == 0, 0, log(df$NHEK.longPolyA.1 + 1))
df$NHEK.longPolyA.norm   = ifelse( df$NHEK.longPolyA   + 1 == 0, 0, log(df$NHEK.longPolyA   + 1))


model.lpa = lm(formula = df$NHEK.longPolyA.norm ~ df$NHEK.longPolyA.norm.1)

ggplot(df, aes(x=log(log(NHEK.longPolyA.1 + 1)),y=log(log(NHEK.longPolyA+ 1)),color=factor(label),size=(3 + label)))+geom_point() + theme_bw()+
  #ggplot(df, aes(x=log(NHEK.longPolyA.norm.1),y=log(NHEK.longPolyA.norm),color=factor(label),size=(3 + label)))+geom_point() + theme_bw()+
  
  #geom_abline(slope= model.lpa$coefficients[2],intercept= model.lpa$coefficients[1]) +
  ggtitle(paste("RNA-seq biological replicates\nlong Poly A w/ psuedocount:(expr = expr + 1)\nCorrelation = ",round(cor( df$NHEK.longPolyA.1 , df$NHEK.longPolyA),digits=5)))
ggsave(file=paste(outdir,"NHEK-longPolyA-compare.pdf"))


#df = df[which(df$NHEK.longNonPolyA.norm.1 * df$NHEK.longNonPolyA.norm > 0), ]

df$NHEK.longNonPolyA.norm.1 =ifelse( df$NHEK.longNonPolyA.1 + 1 == 0, 0, log(log(df$NHEK.longNonPolyA.1 + 1) + 1))
df$NHEK.longNonPolyA.norm   =ifelse( df$NHEK.longNonPolyA   + 1 == 0, 0, log(log(df$NHEK.longNonPolyA   + 1) + 1))

model.lnpa = lm(formula = df$NHEK.longNonPolyA.norm ~ df$NHEK.longNonPolyA.norm.1 )

ggplot(df, aes(x=log(log(NHEK.longNonPolyA.1+1)),y=log(log(NHEK.longNonPolyA+1)),color=factor(label),size=(3 + label)))+geom_point() + theme_bw()+
#ggplot(df, aes(x=NHEK.longNonPolyA.norm.1,y=NHEK.longNonPolyA.norm,color=factor(label),size=(3 + label*1))) + geom_point() + theme_bw()+
  #geom_abline(slope= model.lnpa$coefficients[2],intercept= model.lnpa$coefficients[1])+
 ggtitle(paste("RNA-seq biological replicates\nlong non-Poly A\nCorrelation = ",round(cor( df$NHEK.longNonPolyA.1 , df$NHEK.longNonPolyA),digits=5)))
ggsave(file=paste(outdir,"NHEK-longNonPolyA-compare.pdf"))



plotBothColsZero(df,"NHEK.longNonPolyA.1","NHEK.longNonPolyA","NHEK LongNonPolyA\n found in trial 1 & 2 ?",fileBase=paste(outdir,"NHEK-longNonPolyA"))
plotBothColsZero(df,"NHEK.longPolyA.1","NHEK.longPolyA","NHEK LongPolyA\nfound in trial 1 & 2?",fileBase=paste(outdir,"NHEK-longPolyA"))




}


plotBothColsZero = function(df,col1,col2,title,fileBase){
  c1 = df[[col1]] > 0
  c2 = df[[col2]] > 0
 
  both = c1 & c2
  none = (!c1) & (!c2)
  either = (c1 | c2) & !(c1 & c2)
  df.1 = data.frame(foundIn = c("both","none","either"),geneCount = c(sum(both),sum(none),sum(either)))
  
  ggplot(df.1,aes(x=foundIn,y=geneCount))+geom_histogram(stat="identity")+
    ggtitle(paste(title,"\n",paste(df.1$foundIn,collapse=" "),"\n",paste(df.1$geneCount,collapse=" "),sep=""))
  ggsave(paste(fileBase,"-histogramFound.pdf",sep=""),height=4,width=4)
}




#df.cutoff = df[which(df$NHEK.longPolyA.1 > exp(-4) & df$NHEK.longPolyA > exp(-4)),]
#df.cutoff = df[which(df$NHEK.longPolyA.1 > 0 & df$NHEK.longPolyA > 0),]



test = function(){
 d = data.frame(x = 1:1000/1000 + 100,y=1:1000+runif(1000)*100+10000)
  m = lm(formula=d$y~ d$x )
  ggplot(d,aes(x=x,y=y))+geom_point()+geom_abline(slope=m$coefficients[2],intercept= m$coefficients[1],size=4,color="green")+theme_bw()
  
}

