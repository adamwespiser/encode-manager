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



test = function(){
 d = data.frame(x = 1:1000/1000 + 100,y=1:1000+runif(1000)*100+10000)
  m = lm(formula=d$y~ d$x )
  ggplot(d,aes(x=x,y=y))+geom_point()+geom_abline(slope=m$coefficients[2],intercept= m$coefficients[1],size=4,color="green")+theme_bw()
  
}

