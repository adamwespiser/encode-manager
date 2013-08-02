# Purpose: compare the features of of known lncRNA and the rest of lncRNA
home <- Sys.getenv("HOME")

projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

#boiler plate helpers
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)

#pull in exprLib for some helper functions, and getENSGfromBiomartByRefseq.R for the gene BioMart grabs..
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R")
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R")






#source("./analysis/rnaSeq/compareTopLncsToRest.R")

makeComparisonsLnc <- function(lncDf=df.1,exprColIndex=2:33,transKeyword="annot",
                               fileBase="",outDir="/Users/adam/test/", annotDf, 
                               annotColName,foundColword,verbose=TRUE,
                               injectLncDf = FALSE, processedLncDf){
  #some local helper functions to print debug reports, create outfile names
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  
  
 
  
  
  # prime the df's gene_id
  #lncDf$gene_id_short <- apply(lncDf,1,function(x)as.vector(strsplit(x[c("gene_id")],"\\."))[[1]][1])
  
  # make which genes have been found in annotation

  if (injectLncDf == FALSE){
    lncDf$withinSubset = "false"
    lncDf[which(lncDf$gene_id_short %in% annotDf$gene_id_short),]$withinSubset = "true"
    
    n.subset = length(which(lncDf$withinSubset == "true"))
    n.total = dim(lncDf)[1]
    n.line = paste("subset =",foundColword,"\n",n.subset,foundColword,"lncRNAs out of",n.total,"total lncRNAs\n",sep=" ")
    titleWithBanner = function(x)paste(n.line,x,sep="\n")
    
    
    ddply(annotDf,.(gene_id),function(x)x[1]) -> ensLnc.df
    getLncName <- function(x){y<-ensLnc.df[which(ensLnc.df$gene_id == x ), "lncRnaName"][1];if(is.na(y)){"notFound"}else{y}}
    lncDf[[annotColName]] <- as.vector(sapply(lncDf[["gene_id_short"]],function(x)getLncName(x)))
    #lncDf[is.na(lncDf[[annotColName]]),annotColName] <- "notFound"
    
    print(colnames(lncDf))
    str(lncDf) 
    
  }
  else{lncDf = processedLncDf}
  
  expr.cols <- sort(colnames(lncDf[exprColIndex]))
  expr.uniq.cols <- unique(unlist(lapply(colnames(lncDf[exprColIndex]),function(x)strsplit(x,".long"))))
  expr.cells = expr.uniq.cols[grep("Poly",expr.uniq.cols,invert=TRUE)]
  expr.cols.polyA <-    as.vector(sapply(expr.cells,function(x)paste(x,"longPolyA",sep="."))) 
  expr.cols.nonPolyA <- as.vector(sapply(expr.cells,function(x)paste(x,"longNonPolyA",sep="."))) 
  print(expr.cols.nonPolyA)
  
  
  lpa <- melt(lncDf[c("gene_id","withinSubset",expr.cols.polyA)],measure.vars=sort(expr.cols.polyA),id.vars=c("gene_id","withinSubset"))
  colnames(lpa) <- c("gene_id","withinSubset" ,"expr","longPolyAexpr")
  lpa$cellType <- as.vector(sapply(sapply(as.vector(lpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))
  
  lnpa <- melt(lncDf[c("gene_id","withinSubset",expr.cols.nonPolyA)],measure.vars=sort(expr.cols.nonPolyA),id.vars=c("gene_id","withinSubset"))
  colnames(lnpa) <- c("gene_id","withinSubset", "expr","longNonPolyAexpr")
  lnpa$cellType <- as.vector(sapply(sapply(as.vector(lnpa$expr),function(x){strsplit(x,'\\.')[1]}),function(x)x[1]))
  
  print("a")
  lncDf[["Genc_polyA"]] = "1" 
    pa.sep.df <- merge(merge(lpa,lnpa,by=c("gene_id","cellType","withinSubset")),lncDf[c("Genc_polyA","gene_id","withinSubset")],by=c("gene_id","withinSubset"))
   lpa$seqPullDown <- "longPolyA"
   colnames(lpa) <- c("gene_id", "withinSubset","expr",  "expression" ,"cellType" ,"seqPullDown" )
   lnpa$seqPullDown <- "longNonPolyA"
   colnames(lnpa) <- c("gene_id","withinSubset", "expr",  "expression" ,"cellType" ,"seqPullDown")
   comb <- rbind(lnpa,lpa)
#   
   
   pap.pam.df <- ldply(split(pa.sep.df,pa.sep.df$cellType),getExprComp)
   
  
  p.melt <- melt(pa.sep.df, id.vars=c("gene_id","withinSubset"),measure.var=c("longPolyAexpr","longNonPolyAexpr"))
  p.cell.melt <- melt(pa.sep.df, id.vars=c("gene_id","cellType","withinSubset"),measure.var=c("longPolyAexpr","longNonPolyAexpr"))
  
 
  if (1 == 0){pa.sep.df}
  else {
  printReport("Temp too high...data is melting!!!!")
  id.vec <- c("entropyExpr","threeTimesRest","maxExprType","tissSpec",annotColName,"withinSubset")
  melt.df <- melt(lncDf[c("gene_id",expr.cols,id.vec)], 
                  measure.vars=sort(expr.cols),
                  id.vars=c("gene_id",id.vec))
  
  
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
  
  lncDf$isTissSpec <- "notThreeTimes"
  lncDf[which(lncDf$threeTimesRest != "none"),]$isTissSpec <- "threeTimes"
  tmp <- lncDf
  tmp$isTissSpec <- "combinedDistro"
  tissSpec.df <- rbind(lncDf,tmp)
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
  
  maxExprTable <- as.data.frame(table(lncDf$maxExprType))
  threeTimeTable <- as.data.frame(table(subset(melt.df, threeTimesRest == variable)[c("threeTimesRest")]))
  
##############################  ##############################  ##############################  ##############################
  lncDf$annotName <-with(lncDf,get(annotColName))
  melt.df$annotName <- with(melt.df,get(annotColName))
  measure <- c("sumExpr","varianceExpr","averageExpr","minExpr","maxExpr","nTimesRest","tissSpec")
  stat.melt <- melt(lncDf[c("gene_id",measure,"withinSubset")],measure.vars=measure,id.vars=c("gene_id","withinSubset"))                    
  derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  der.cols <- c("No_Exons","LncRNA_tx_size","phastcons_primate_intron_score","phastcons_primate_transcript_score","phastcons_primate_promoter_score","phastcons_mammal_transcript_score","phastcons_mammal_intron_score","phastcons_mammal_promoter_score","phastcons_vertebrate_transcript_score","phastcons_vertebrate_intron_score","phastcons_vertebrate_promoter_score")
  derr.df <- ddply(derrien.df[c("LncRNA_GeneId",der.cols)],.(LncRNA_GeneId),numcolwise(mean))
  combinedDerr.df <- merge(lncDf[c("gene_id","withinSubset")],derr.df,by.x="gene_id",by.y="LncRNA_GeneId")
  derr.melt <- melt(combinedDerr.df,id.vars=c("gene_id","withinSubset"))
  comDer.df <- merge(combinedDerr.df,lncDf[c("gene_id","lncRnaName")],by="gene_id")
  
  
  allPlotRun <- 1
  if (1 == 1){
    #poly A pulldown plots
    ggplot(p.melt,aes(x=log(value),fill=variable))+geom_density(alpha=I(0.4))+facet_wrap(~withinSubset)+theme_bw()+
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("log(RPKM of individual sample)")
    ggsave(file=makeOutFile("lncRNA-RPKM-pulldown-density-plot.pdf"),height=8,width=8)
    
    ggplot(p.melt,aes(x=log(value),color=variable,fill=variable))+geom_freqpoly(binwidth=1)+facet_wrap(~withinSubset,scale="free_y")+theme_bw()+
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("log(RPKM of individual sample)")
    ggsave(file=makeOutFile("lncRNA-RPKM-pulldown-freqpoly-plot.pdf"),height=8,width=8)
    
    ggplot(p.cell.melt,aes(x=log(value),fill=variable))+geom_density(alpha=I(0.4))+facet_grid(cellType~withinSubset)+theme_bw()+
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("log(RPKM of individual sample)")
    ggsave(file=makeOutFile("lncRNA-RPKM-cellType-subset-density-plot.pdf"),height=14,width=8)
    
    ggplot(p.cell.melt,aes(x=log(value),color=variable,fill=variable))+geom_freqpoly(binwidth=1)+facet_grid(cellType~withinSubset)+theme_bw()+
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("log(RPKM of individual sample)")
    ggsave(file=makeOutFile("lncRNA-RPKM-cellType-subset-freqpoly-plot.pdf"),height=14,width=8)
    
    # what is the average expression of known lnc's versus the rest? 
    ggplot(melt.df,aes(x=log(value),fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("log(RPKM of individual sample)")
    ggsave(file=makeOutFile("lncRNA-RPKM-density-plot.pdf"),height=8,width=8)
    
    #what is the JSD values between found in disease and otherwise?
    ggplot(lncDf,aes(x=tissSpec,fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("JSD for each gene")
    ggsave(file=makeOutFile("lncRNA-JSD-density-plot.pdf"),height=8,width=8)
    
    ggplot(lncDf,aes(x=tissSpec,fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+ facet_wrap(~withinSubset)
      ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
      xlab("JSD for each gene")
    ggsave(file=makeOutFile("lncRNA-JSD-density-plot.pdf"),height=8,width=8)
    
    
    ggplot(lncDf,aes(x=log(tissSpec),y=log(nTimesRest)))+geom_point()+ theme_bw()+facet_wrap(~withinSubset)+
      ggtitle(titleWithBanner("tissue spec. vs. nTimesRest expression measure\nFacet by foundInDisease?"))+
      xlab("JSD for each Gene")+ylab("nTimesRest")
    ggsave(file=makeOutFile("lncRNA-nTimesRest-JSD-scatter-plot.pdf"),height=8,width=8)
    
    
    if (injectLncDf == FALSE){
    ggplot(melt.df[which(melt.df$withinSubset == "true"),],aes(x=variable,y=log2(value),group=gene_id,color=tissSpec)) + 
      geom_line()+coord_polar()+
      theme_bw()+
      ggtitle(titleWithBanner("star plot of total expression over all samples"))
    ggsave(file=makeOutFile("subsetOnly-polar.pdf"),height=16,width=16)
    
    ###
    #
    #             COORD POLAR
    #
    ggplot(melt.df,aes(x=variable,y=log2(value),group=gene_id,color=withinSubset,alpha=0.6)) + geom_line()+
      geom_line(data=melt.df[which(melt.df$withinSubset == "true"),],aes(x=variable,y=log2(value),group=gene_id,color="true"))+
      
      coord_polar()+theme_bw()+
      scale_colour_manual(values = c("false"="lightgrey","true"="blue","red"="blue"))  +
      ggtitle(titleWithBanner("star plot of total expression over all samples"))
    ggsave(file=makeOutFile("allEntries-ColorBySubset-polar.pdf"),height=16,width=16)
    
    }
    #  id.vec <- c("entropyExpr","threeTimesRest","maxExprType","tissSpec",annotColName,foundColword)
    #  melt.df <- melt(lncDf[c("gene_id",expr.cols,id.vec)], 
    #                  measure.vars=sort(expr.cols),
    #                  id.vars=c("gene_id",id.vec))
    
    
    ggplot(stat.melt,aes(x=log(value),fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
      facet_wrap(~variable,scale="free") +
      ggtitle(titleWithBanner("measures accross all samples(cell-type/pulldown)"))+
      xlab("value of facet label")+
      scale_colour_manual(values = c("true"="green","false"="red"))
    ggsave(file=makeOutFile("lncRNA-allcols-density-plot.pdf"),height=8,width=8)
    
    ggplot(stat.melt,aes(x=log(value),fill=withinSubset))+geom_freqpoly() + theme_bw()+
      facet_grid(withinSubset~variable,scale="free_y") +
      ggtitle(titleWithBanner("measures accross all samples(cell-type/pulldown)"))+
      xlab("value of facet label")+
      scale_colour_manual(values = c("true"="green","false"="red"))
    ggsave(file=makeOutFile("lncRNA-allcols-freqpoly-plot.pdf"),height=6,width=13)
    
    
    ggplot(stat.melt,aes(x=withinSubset,y=log(value),fill=withinSubset))+geom_boxplot() + theme_bw()+
      facet_wrap(~variable,scale="free") +
      ggtitle(titleWithBanner("measures accross all samples(cell-type/pulldown)"))+
      xlab("value of facet label")
    #+scale_colour_manual(values = c("true"="green","false"="red"))
    ggsave(file=makeOutFile("lncRNA-allcols-box-plot.pdf"),height=8,width=8)
    
    
    
    
    ggplot(derr.melt,aes(x=log(value),fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
      facet_wrap(~variable,scale="free") +
      ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
      xlab("value of facet label")+
      scale_colour_manual(values = c("true"="green","false"="red"))
    ggsave(file=makeOutFile("lncRNA-derrien-allcols-density-plot.pdf"),height=8,width=8)
    
    ggplot(derr.melt,aes(x=log(value),fill=withinSubset))+geom_freqpoly() + theme_bw()+
      facet_grid(withinSubset~variable,scale="free_y") +
      ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
      xlab("value of facet label")+
      scale_colour_manual(values = c("true"="green","false"="red"))
    ggsave(file=makeOutFile("lncRNA-derrien-allcols-freqpoly-plot.pdf"),height=6,width=16)
    
    ggplot(derr.melt,aes(x=withinSubset,y=log(value),fill=withinSubset))+geom_boxplot() + theme_bw()+
      facet_wrap(~variable,scale="free") +
      ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
      xlab("value of facet label")
    #+scale_colour_manual(values = c("true"="green","false"="red"))
    ggsave(file=makeOutFile("lncRNA-derrien-allcols-box-plot.pdf"),height=8,width=8)
    
    #lets run some wilcox tests...
    if (injectLncDf == FALSE){
      
    print("starting statistical comparisons (pairwise)")
    allPairs.df <- merge(combinedDerr.df,lncDf[c("gene_id",measure)],by="gene_id")
    within.df  <- allPairs.df[allPairs.df$withinSubset == "true",]
    without.df <- allPairs.df[allPairs.df$withinSubset == "false",]
    wilcoxFun <- function(x)wilcox.test(within.df[,c(x)],without.df[,c(x)])$p.value
    ksFun <- function(x)ks.test(within.df[,c(x)],without.df[,c(x)])$p.value
    tFun <- function(x)t.test(within.df[,c(x)],without.df[,c(x)])$p.value
    #partDF <- function(x)data.frame(col=x,wilcox.pval=wilcoxFun(x),ks.pval=ksFun(x),ttest.pval=tFun(x))
    partDF <- function(x)data.frame(col=x,wilcox.pval=wilcoxFun(x))
    
    pvals.df <- ldply(as.list(c(measure,der.cols)),partDF)
    
    ggplot(pvals.df, aes(x=col,y=wilcox.pval))+geom_histogram()+geom_abline(slope=0,intercept=0.05,color="red")+coord_flip()+
      ylab("pValue from wilcox test")+
      xlab("transcript feature")+
      ggtitle(titleWithBanner("pvalues from Derrien 2012 datasheet and expression data"))
    ggsave(file=makeOutFile("lncRNA-derrien-allcols-pvalues.pdf"),height=8,width=8)
    
    ggplot(pvals.df, aes(x=col,y=log(wilcox.pval)))+geom_histogram()+geom_abline(slope=0,intercept=0.05,color="red")+coord_flip()+
      ylab("log(pValue from wilcox test)")+
      xlab("transcript feature")+
      ggtitle(titleWithBanner("wilcox pvalues from Derrien 2012 datasheet and expression data"))
    ggsave(file=makeOutFile("lncRNA-derrien-allcols-LOGpvalues.pdf"),height=8,width=8)
    }
  } 
  
dispersion.df <- merge(lncDf,combinedDerr.df[c("gene_id",der.cols)],by="gene_id")
  
  
  disp.cols <- c(der.cols,measure)
  
  
makeDispersionPlot <- function(plot.df,colName){
 outSapplyLocal <- paste("lncRNA-",colName,"-dispersion-plot-subSetLabelled.pdf",sep="")

 x.high = max(log(plot.df[[colName]])+2)
 x.low =  min(log(plot.df[[colName]])+4)
 
 plot.df[["hjust"]] <- runif(dim(plot.df)[1]) -0.5
 #    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
 #    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
 colFunStr <- paste("log(",colName,")",sep="")
 found.df <- plot.df[plot.df$lncRnaName != "notFound",]
 found.df[which(is.finite(log(found.df[[colName]] + 1)) != 0),]
 
 notFound.df <- plot.df[plot.df$lncRnaName == "notFound",]
 print("intial munge")
 labels.df <- getHjustForRange(log(as.numeric(found.df[[colName]])),found.df[["lncRnaName"]])
 str(labels.df)
 
 n.count <- paste(dim(labels.df)[1],"in subset out of ", dim(plot.df)[1],"total lncRNAs\n")
 title <- paste("jitter plot of ",colName," data\n","subsetted lncs w/ label\n",n.count,sep="")
 print("going for ggplots")
 ggplot(notFound.df)+
    geom_jitter(alpha=I(0.2),width=1,aes_string(x=colFunStr,y="0",color="withinSubset"))+
    coord_flip()+theme_bw()+
    layer(data=labels.df,mapping=aes(x=x,y=y,label=label),geom="text",size=3)+
    annotate("text",x=x.high,y=-0.5,label="(high)")+
    annotate("text",x=x.low,y=-0.5,label="(low)")+
    scale_color_manual(values = c("true"="black","false"="red"))+
    ggtitle(titleWithBanner(title))+
    xlab("x =  random uniform[-0.5,0,5]")
  ggsave(file=makeOutFile(outSapplyLocal),height=8,width=8)
}

  if (dim(dispersion.df[dispersion.df$lncRnaName != "notFound",])[1] < 100){
sapply(disp.cols,function(x)makeDispersionPlot(colName=x,plot.df=dispersion.df ))
  }  

  #dispersion.df
} # extra "else"

}

calcOffset <-function(x){if(runif(1) < 0.5){ seq(0,1,by=1/length(x$vValue))[seq_along(x$vValue)]}else{ 1-seq(0,1000,by=1/length(x$vValue))[seq_along(x$vValue)]}}

getHjustForRange <- function(values,names){
  notInf <- which(is.finite(values)) 
  values <- values[notInf]
  names <- names[notInf]
  
  col.min <-min(values)
   col.max <- max(values)
   #print(c(col.min,col.max))
   range<- seq(to=col.max,from=col.min,by=abs(col.max-col.min)/15)
   ddply(data.frame(vRange= cut(values,breaks=range),vValue=values,labels=names),
         .(vRange),
         function(x)
           data.frame(label=x$labels,
                      x=x$vValue,
                      y=((calcOffset(x)-0.5)*0.9)))
                      #x=seq(0,1000,by=0.11)[seq_along(x$vValue)]))
}


runComps <- function(x){
  #grab the (gene = top transcript over sample) for lncRNA data. 
  combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
  lnc.out.dir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots_eachSample"
  lnc.fileBase = "lncRNA-gene=MaxTransSample"
  lnc.transKeyword = "lncRNA genes as max trans in sample"
  lnc.colIndex = 2:33
  df.1=readInTable(lnc.in.file)
  df.1$gene_id_short <- sapply(df.1$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
  
  # Get lncRNAs that are 1) known to be involved in disease 2) Identified through literature search as "functional verified and mechanism determined"
  lncDisease.df <- getlncDiseaseDf()# colnames = c("lnRnaName,refseq,gene_id)
  lncDisease.df$gene_id_short <- lncDisease.df$gene_id
  lncDisease.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncDiseaseDb"
  print("start lnc Disease")
  makeComparisonsLnc(lncDf=df.1,exprColIndex=2:33,foundColword="foundInDisease",fileBase="lncRnaDisease-",
                     outDir=lncDisease.outdir,annotDf=lncDisease.df,annotColName="lncRnaName")
  
  
  lncdList.df <- getlncdListDf()# colnames = c("lnRnaName,refseq,gene_id)
  lncdList.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaDbPlusLit"
  print("start lnc List")
  
  lncdList.df[["gene_id_short"]] <- sapply(lncdList.df$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  lncdList.df[["gene_id"]] <-lncdList.df[["gene_id_short"]]
  makeComparisonsLnc(lncDf=df.1,exprColIndex=2:33,foundColword="foundLiterature",fileBase="lncRnaLiterature-",
                     outDir=lncdList.outdir,annotDf=lncdList.df[lncdList.df$foundInDerrien == "yes",],annotColName="lncRnaName")
  

  func.df <- getEnslist()
  #func1.df <- readInTable(getFullPath("data/lncRnaExpr_ml.tab"))
  
  func.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaFunctional"
  print("start functional List")
  func.df[["gene_id_short"]] <-  func.df[["ensembl_gene_id"]] 
  func.df[["gene_id"]] <-func.df[["gene_id_short"]]
  func.df[["lncRnaName"]] <- func.df[["external_gene_id"]]
  makeComparisonsLnc(lncDf=df.1,exprColIndex=2:33,foundColword="functionalLncRNA",fileBase="func-",
                     outDir=func.outdir,annotDf=unique(func.df[c( "lncRnaName","gene_id","gene_id_short")]) ,annotColName="lncRnaName")
  
  
  lncMouseOrthos.df <- getOrthosDf() # "lncRnaName"     "gene_id"        "foundInDerrien"
  lncMouseOrthos.df$gene_id <- as.character(lncMouseOrthos.df$gene_id)
  lncMouseOrthos.df$gene_id_short <- sapply(lncMouseOrthos.df$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  lncMouseOrthos.df$gene_id <-  lncMouseOrthos.df$gene_id_short
  lncMouseOrthos.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncOrthoToMouse"
  print("start mouse orthos")
  makeComparisonsLnc(lncDf=df.1,exprColIndex=2:33,foundColword="foundInMouse",fileBase="mouse-ortho-",
                     outDir=lncMouseOrthos.outdir,annotDf=lncMouseOrthos.df,annotColName="lncRnaName")
  
  
  
  
  antiSense.df <- getDerriensAntiSenseDf()
  antiSense.df$gene_id <- as.character(antiSense.df$gene_id)
  antiSense.df$gene_id_short <- sapply(antiSense.df$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  lncMouseOrthos.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/antiSenseByDerrien"
  print("start derrien AS")
  makeComparisonsLnc(lncDf=df.1,exprColIndex=2:33,foundColword="antisenseGene",fileBase="as-gene-",
                     outDir=lncMouseOrthos.outdir,annotDf=antiSense.df,annotColName="lncRnaName")
  
  titles <- c("anti-sense", "lncMouseOrthologs","lncRnaDbPlusLit","lncDiseaseDb")
  counts <- c(dim(antiSense.df)[1],dim(lncMouseOrthos.df)[1],dim(lncdList.df)[1],dim(lncDisease.df)[1])
  
  write.csv(file="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/summary.csv",
            data.frame(title=titles,count=counts))
  
}




getAnnotLncDf <- function(lncDf,annotDf,exprCol,annotColName){

  lncDf$withinSubset = "false"
  lncDf[which(lncDf$gene_id_short %in% annotDf$gene_id_short),]$withinSubset = "true"
  
  #apply annotDf information to lncDf
  ddply(annotDf,.(gene_id),function(x)x[1,]) -> ensLnc.df
  getLncName <- function(x){ensLnc.df[which(ensLnc.df$gene_id == x ), "lncRnaName"][1]}
  lncDf[c(annotColName)] <- apply(lncDf,1,function(x)getLncName(x[c("gene_id_short")]))
  lncDf[is.na(lncDf[c(annotColName)]),annotColName] <- "notFound"
  lncDf
}

withinFuntionalCompareBiotypes = function(){
  
  
lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
lncFound.df = getEnslist()
lncExpr.df = readInTable(lnc.in.file)
lncExpr.df[["gene_id_short"]] =  sapply(as.character(lncExpr.df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
lncFound.df = unique(data.frame(ensembl_gene_id=ens$ensembl_gene_id,external_gene_id=ens$external_gene_id,gene_biotype=ens$gene_biotype))
lncFound.df = lncFound.df[which(lncFound.df$ensembl_gene_id %in% lncExpr.df[["gene_id_short"]]), ]
lncFoundThree.df = lncFound.df[which(lncFound.df$gene_biotype %in% c("antisense","lincRNA","processed_transcript")),]

lncExprFound.df = lncExpr.df[which(lncExpr.df$gene_id_short %in% lncFoundThree.df$ensembl_gene_id),]
lncExprFound.df = merge(lncExprFound.df,lncFound.df,by.x="gene_id_short",by.y="ensembl_gene_id")

exprCols = colnames(lncExprFound.df)[3:34]

lncExprFound.df$withinSubset = lncExprFound.df$gene_biotype
annotColName = "external_gene_id"
lncExprFound.df$lncRnaName = lncExprFound.df[[annotColName]]


lncExprFound.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncFunctional-BioTypes"
print("start biotype")
makeComparisonsLnc(lncDf=lncExpr.df,exprColIndex=exprCols,foundColword="biotypeOfFunctional",fileBase="func-biotype-",
                   outDir=lncExprFound.outdir,annotDf=lncExprFound.df,annotColName="external_gene_id",
                   injectLncDf = TRUE, processedLncDf=lncExprFound.df)


lncExpr.df$external_gene_id = lncExpr.df$gene_id_short 
lncExpr.df$withinSubset = "rest"
lncExpr.df$lncRnaName = lncExpr.df[[annotColName]]
lncExpr.df$gene_biotype = "nonFunctional"
lncExpr.df$lncRnaName = lncExpr.df$gene_id_short
lncExprMod.df = lncExpr.df[which(!lncExpr.df$gene_id_short %in% lncExprFound.df$gene_id_short),]
lncExprMod.df = lncExprMod.df[colnames(lncExprFound.df)]
lncExprFound.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncFunctional-BioTypesVsRest"
print("start biotype")
makeComparisonsLnc(lncDf=lncExpr.df,exprColIndex=exprCols,foundColword="biotypeOfFunctional",fileBase="func-biotype-",
                   outDir=lncExprFound.outdir,annotDf=lncExprFound.df,annotColName="external_gene_id",
                   injectLncDf = TRUE, processedLncDf=rbind(lncExprMod.df,lncExprFound.df))


geneTxType.df = getBioTypesForDf()

geneTxType.df = unique(derrien.df[c("LncRNA_GeneId","tx_biotype")])
geneTxType.df[["gene_id_short"]] =  sapply(as.character(geneTxType.df$LncRNA_GeneId),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
geneTxTypeFunc.df = geneTxType.df[which(geneTxType.df$LncRNA_GeneId %in% lncExprFound.df$gene_id),]
lncExprBmBioTypes.df = unique(merge(lncExpr.df,geneTxType.df[c("LncRNA_GeneId","bm_biotype")],by.x="gene_id",by.y="LncRNA_GeneId"))



#antisense
lncExprBmBioTypes.df$withinSubset = "unlabel-AS"
lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$gene_id_short %in% lncFound.df$ensembl_gene_id),"withinSubset"] = "functional-AS"
lncExprFound.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncFunctional-BioTypes-AS"
makeDir(lncExprFound.outdir)
makeComparisonsLnc(lncDf=lncExpr.df,exprColIndex=exprCols,foundColword="antisense",fileBase="func-biotype-",
                   outDir=lncExprFound.outdir,annotDf=lncExprFound.df,annotColName="external_gene_id",
                   injectLncDf = TRUE, processedLncDf= lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$bm_biotype == "antisense"), ])

#processed_transcript
lncExprBmBioTypes.df$withinSubset = "unlabel-procTrans"
lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$gene_id_short %in% lncFound.df$ensembl_gene_id),"withinSubset"] = "functional-procTrans"
lncExprFound.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncFunctional-BioTypes-PT"
makeDir(lncExprFound.outdir)
makeComparisonsLnc(lncDf=lncExpr.df,exprColIndex=exprCols,foundColword="processed_transcript",fileBase="func-processed-transcript",
                   outDir=lncExprFound.outdir,annotDf=lncExprFound.df,annotColName="external_gene_id",
                   injectLncDf = TRUE, processedLncDf= lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$bm_biotype == "processed_transcript"), ])



#lincRNA
lncExprBmBioTypes.df$withinSubset = "unlabel-lincRNA"
lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$gene_id_short %in% lncFound.df$ensembl_gene_id),"withinSubset"] = "functional-lincRNA"
lncExprFound.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncFunctional-BioTypes-lincRNA"
makeDir(lncExprFound.outdir)
makeComparisonsLnc(lncDf=lncExpr.df,exprColIndex=exprCols,foundColword="lincRNA",fileBase="func-lincRNA-",
                   outDir=lncExprFound.outdir,annotDf=lncExprFound.df,annotColName="external_gene_id",
                   injectLncDf = TRUE, processedLncDf= lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$bm_biotype == "lincRNA"), ])



threeTypesTotal.df <- lncExprBmBioTypes.df[which(lncExprBmBioTypes.df$bm_biotype %in% c("lincRNA","processed_transcript","antisense")),]
threeTypesTotal.df$functional <- "unLabel-"
threeTypesTotal.df[which(threeTypesTotal.df$gene_id_short %in% lncFound.df$ensembl_gene_id),"functional"] <- "func-"
threeTypesTotal.df$withinSubset = paste(threeTypesTotal.df$bm_biotype,threeTypesTotal.df$functional,sep="")
lncExprFound.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncFunctional-BioTypes-fullCompare"
makeDir(lncExprFound.outdir)
makeComparisonsLnc(lncDf=lncExpr.df,exprColIndex=exprCols,foundColword="funcBioType",fileBase="func-biotype-",
                   outDir=lncExprFound.outdir,annotDf=lncExprFound.df,annotColName="external_gene_id",
                   injectLncDf = TRUE, processedLncDf= threeTypesTotal.df)


}



getBioTypesForDf <- function(df=geneTxType.df,file=getFullPath("data/lnc_biotype.tab")){
  if (!file.exists(file)){
    numberGroups = 50
    df$group = cut(seq_along(df$LncRNA_GeneId),breaks=50,labels = 1:numberGroups)
  geneBiotypes = sapply(1:numberGroups,function(x){sapply(geneTxType.df[which(df$group == x),"gene_id_short"],
                                                          function(x){
                                                            x = getBM(attributes=c("gene_biotype"),filters=c("ensembl_gene_id"),values=list(x),mart=mart)
                                                            if (is.logical(x)){
                                                              return("notfound")
                                                            }
                                                            return(x)
                                                          }   )})
    df$bm_biotype = unlist(geneBiotypes)
  #lnc.biotype.out = 
  exportAsTable(file=file,df=geneTxType.df)
  } else{
     d = readInTable(file=file) 
  } 
  
}


pcaAnalysis <- function(lncDf,exprCols,normalFun=FALSE,outDir,fileBase,foundColword){
  
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner = function(x)paste(paste("lncRNA subset = ",foundColword),x,sep="\n")
  expr.cols <- exprCols
  
  
  lnc.pca             <- princomp(lncDf[exprCols],cor=FALSE)
  #lnc.pca.cov         <- princomp(covmat= cov(lncDf[exprCols]))
  lnc.pca.cor         <- princomp(lncDf[exprCols],cor=TRUE)
  
  #max(abs(dist(lncDf[expr.cols]),dist(cmdscale(dist(lncDf[expr.cols]), k= 5))))
  
  data.frame(pca.std.dev=lnc.pca$sdev^2,
                    cor.std.dev=lnc.pca.cor$sdev^2,
                   # cov.std.dev=lnc.pca.cov$sdev^2,
                    component=1:32) -> comp.var
  comp.var.melt <- melt(comp.var, id.var="component")  
  
  ggplot(comp.var.melt,aes(x=component,y=value))+geom_line()+theme_bw()+facet_wrap(~variable)+
    ggtitle("PCA lncRNA expressionmatrix\nComponent variation over all components")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-egeinvals.pdf"),height=8,width=8)
  
  ggplot(data.frame(cumVar=c(0,cumsum(comp.var$pca.std.dev^2)/sum(comp.var$pca.std.dev^2)),component=0:32),aes(x=component,y=cumVar))+geom_line()+geom_point()+theme_bw()+
    ggtitle("PCA lncRNA expressionmatrix\nComponent variation over all components")+xlab("cummulative component")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-egeinvals-ecdf.pdf"),height=8,width=8)
  
  ggplot(data.frame(cumVar=c(0,cumsum(comp.var$cor.std.dev^2)/sum(comp.var$cor.std.dev^2)),component=0:32),aes(x=component,y=cumVar))+geom_line()+geom_point()+theme_bw()+
    ggtitle("PCA lncRNA CORRELATION expressionmatrix\nCumulative sum of component")+xlab("cummulative component")
  ggsave(file=makeOutFile("lncRNA-pca-cor-exprVectors-egeinvals-ecdf.pdf"),height=8,width=8)
  
  
  lnc.pca.factors <-  as.data.frame(lnc.pca$scores)
  lnc.pca.factors[["gene_id"]] <- lncDf[["gene_id"]]
  lnc.pca.factors[["withinSubset"]] <-lncDf[["withinSubset"]]
  lnc.pca.factors[["lncRnaName"]] <-lncDf[["lncRnaName"]]
  
  lnc.pca.factors[["Comp.1"]] <-lnc.pca.factors[["Comp.1"]] - min(lnc.pca.factors[["Comp.1"]]) + 1
  lnc.pca.factors[["Comp.2"]] <-lnc.pca.factors[["Comp.2"]] - min(lnc.pca.factors[["Comp.2"]]) + 1
  lnc.pca.factors[["Comp.3"]] <-lnc.pca.factors[["Comp.3"]] - min(lnc.pca.factors[["Comp.3"]]) + 1
  lnc.pca.factors[["Comp.4"]] <-lnc.pca.factors[["Comp.4"]] - min(lnc.pca.factors[["Comp.4"]]) + 1
  
  lnc.pca.factors.cor <-  as.data.frame(lnc.pca.cor$scores)
  lnc.pca.factors.cor[["gene_id"]] <- lncDf[["gene_id"]]
  lnc.pca.factors.cor[["withinSubset"]] <- lncDf[["withinSubset"]]
  lnc.pca.factors.cor[["lncRnaName"]] <-lncDf[["lncRnaName"]]
  
  lnc.pca.factors.cor[["Comp.1"]] <-lnc.pca.factors.cor[["Comp.1"]] - min(lnc.pca.factors.cor[["Comp.1"]]) + 1
  lnc.pca.factors.cor[["Comp.2"]] <-lnc.pca.factors.cor[["Comp.2"]] - min(lnc.pca.factors.cor[["Comp.2"]]) + 1
  lnc.pca.factors.cor[["Comp.3"]] <-lnc.pca.factors.cor[["Comp.3"]] - min(lnc.pca.factors.cor[["Comp.3"]]) + 1
  lnc.pca.factors.cor[["Comp.4"]] <-lnc.pca.factors.cor[["Comp.4"]] - min(lnc.pca.factors.cor[["Comp.4"]]) + 1
   
  PCA.title <- titleWithBanner( "PCA of [expression vectors] foreach lnc ")
  ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=3)+
    ggtitle(PCA.title) +
    xlab("log(First Component)") + ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors.pdf"),height=8,width=8)
  
  
  
  PCA.title.cor <- titleWithBanner( "PCA of CORRELATION [expression vectors] foreach lnc ")
  ggplot(lnc.pca.factors.cor,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=3)+
    ggtitle(PCA.title.cor) +
    xlab("log(First Component)") + ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-cor.pdf"),height=8,width=8)
 
  ggplot(lnc.pca.factors.cor,aes(x=log(Comp.1),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.3)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.3),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=3)+
    ggtitle(PCA.title.cor) +
    xlab("log(First Component)") + ylab("log(Third Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-cor-1-3.pdf"),height=8,width=8)  
  
  ggplot(lnc.pca.factors.cor,aes(x=log(Comp.2),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.2),y=log(Comp.3)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.2),y=log(Comp.3),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(PCA.title.cor) +
    xlab("log(Second Component)") + ylab("log(Third Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-cor-2-3.pdf"),height=8,width=8)
  
  
  lnc.sub <- which(lnc.pca.factors$withinSubset == "true")
  lnc.min <- lnc.pca.factors[which(lnc.pca.factors[["Comp.1"]] != min(lnc.pca.factors[["Comp.1"]])
                                   & lnc.pca.factors[["Comp.2"]] != min(lnc.pca.factors[["Comp.2"]])
                                   & lnc.pca.factors[["Comp.3"]] != min(lnc.pca.factors[["Comp.3"]]) ),]
  
  ggplot(lnc.min,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title,"\nwithout outlier"))+
    xlab("log(First Component)")+
    ylab("log(Second Component)")
    ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min.pdf"),height=8,width=8)
  
  ggplot(lnc.min,aes(x=log(Comp.1),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title,"\nwithout outlier"))+
    xlab("log(First Component)")+
    ylab("log(Third Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min-comp1-3.pdf"),height=8,width=8)
  
  ggplot(lnc.min,aes(x=log(Comp.2),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title,"\nwithout outlier"))+
    xlab("log(Second Component)")+
    ylab("log(Third Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min-comp2-3.pdf"),height=8,width=8)
  
  lnc.sub.cor <- which(lnc.pca.factors.cor$withinSubset == "true")
  lnc.min.cor <- lnc.pca.factors.cor[which(lnc.pca.factors.cor[["Comp.1"]] != min(lnc.pca.factors.cor[["Comp.1"]])),]
  lnc.min <- lnc.pca.factors.cor[which(lnc.pca.factors.cor[["Comp.1"]] != min(lnc.pca.factors.cor[["Comp.1"]])
                                   & lnc.pca.factors.cor[["Comp.2"]] != min(lnc.pca.factors.cor[["Comp.2"]])
                                   & lnc.pca.factors.cor[["Comp.3"]] != min(lnc.pca.factors.cor[["Comp.3"]]) ),]
  
  ggplot(lnc.min.cor,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title.cor,"\nwithout outlier CORRELATION"))+
    xlab("log(First Component)")+
    ylab("log(Second Component)")
    ggsave(file=makeOutFile("lncRNA-pca-cor-exprVectors-remove-min.pdf"),height=8,width=8)
  
  ggplot(lnc.min.cor,aes(x=log(Comp.1),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title.cor,"\nwithout outlier CORRELATION"))+
    xlab("log(First Component)")+
    ylab("log(Third Component)")
  ggsave(file=makeOutFile("lncRNA-pca-cor-exprVectors-remove-min-com1-3.pdf"),height=8,width=8)
  
  ggplot(lnc.min.cor,aes(x=log(Comp.2),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title.cor,"\nwithout outlier CORRELATION"))+
    xlab("log(Second Component)")+
    ylab("log(Third Component)")
  ggsave(file=makeOutFile("lncRNA-pca-cor-exprVectors-remove-min-comp2-3.pdf"),height=8,width=8)
  
  
  # lnc.min.cor lnc.min    lnc.pca.factors.cor lnc.pca.factors
  
  runStatHexBin <- function(df,title,fileBase){
    ggplot(df,aes(x=log(Comp.1),y=log(Comp.2))) + stat_binhex()+theme_bw()+ggtitle(paste(title))
    ggsave(file=makeOutFile(paste(fileBase,"-1-2.pdf",sep="")),height=8,width=8)
    
    ggplot(df,aes(x=log(Comp.1),y=log(Comp.3))) + stat_binhex()+theme_bw()+ggtitle(title)
    ggsave(file=makeOutFile(paste(fileBase,"-1-3.pdf",sep="")),height=8,width=8)
    
    ggplot(df,aes(x=log(Comp.2),y=log(Comp.3))) + stat_binhex()+theme_bw()+ggtitle(title)
    ggsave(file=makeOutFile(paste(fileBase,"-2-3.pdf",sep="")),height=8,width=8)
    }
  
  runStatHexBin(df=lnc.min.cor,title="lnc correlation PCA\nmin removed",fileBase="lncRNA-cor-nomin-binhex")
  runStatHexBin(df=lnc.min,title="lnc PCA\nmin removed",fileBase="lncRNA-nomin-binhex")
  runStatHexBin(df=lnc.pca.factors.cor,title="lnc correlation PCA",fileBase="lncRNA-cor-binhex")
  runStatHexBin(df=lnc.pca.factors,title="lnc PCA",fileBase="lncRNA-binhex")
  
}


runPCA <- function(x){
  #grab the (gene = top transcript over sample) for lncRNA data. 
  combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
  lnc.out.dir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots_eachSample"
  lnc.fileBase = "lncRNA-gene=MaxTransSample"
  lnc.transKeyword = "lncRNA genes as max trans in sample"
  lnc.colIndex = 2:33
  df.1=readInTable(lnc.in.file)
  df.1$gene_id_short <- sapply(df.1$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
  
  # Get lncRNAs that are 1) known to be involved in disease 2) Identified through literature search as "functional verified and mechanism determined"
  lncDisease.df <- getlncDiseaseDf()# colnames = c("lnRnaName,refseq,gene_id)
  lncDisease.df$gene_id_short <- lncDisease.df$gene_id
  lncDisease.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncDiseaseDb"
  print("start lnc Disease")
  lncDisease.df <- getAnnotLncDf(lncDf=df.1,annotDf=lncDisease.df,exprCol=2:33,annotColName="lncRnaName")
  pcaAnalysis(lncDf=lncDisease.df,exprCols=lnc.colIndex,normalFun=FALSE,foundColword="foundInDisease",outDir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncDiseaseDb",fileBase="PCA-")
  

  
  lncdList.df <- getlncdListDf()# colnames = c("lnRnaName,refseq,gene_id)
  lncdList.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaDbPlusLit"
  print("start lnc List")
  
  lncdList.df[["gene_id_short"]] <- sapply(lncdList.df$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  lncdList.df[["gene_id"]] <-lncdList.df[["gene_id_short"]]
  lncList.df <- getAnnotLncDf(lncDf=df.1,annotDf=lncdList.df,exprCol=2:33,annotColName="lncRnaName")
  pcaAnalysis(lncDf=lncList.df,exprCols=lnc.colIndex,normalFun=FALSE,foundColword="foundInLiterature",
              outDir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaDbPlusLit",fileBase="PCA-")
  
  
  func.df <- getEnslist()
  #func1.df <- readInTable(getFullPath("data/lncRnaExpr_ml.tab"))
  
  func.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaFunctional"
  print("start functional List")
  func.df[["gene_id_short"]] <-  func.df[["ensembl_gene_id"]] 
  func.df[["gene_id"]] <-func.df[["gene_id_short"]]
  func.df[["lncRnaName"]] <- func.df[["external_gene_id"]]
  f.df <- getAnnotLncDf(lncDf=df.1,annotDf=func.df,exprCol=2:33,annotColName="lncRnaName")
  pcaAnalysis(lncDf=f.df,exprCols=2:33,foundColword="functional-LncRNA",fileBase="PCA-",
                     outDir=func.outdir)
  

}


