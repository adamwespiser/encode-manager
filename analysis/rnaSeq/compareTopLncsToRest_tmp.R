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




makeShortGeneId <- function(df){
  
  df$
  
  df
}

#source("./analysis/rnaSeq/compareTopLncsToRest.R")

makeComparisonsLnc <- function(lncDf=df.1,exprColIndex=2:33,transKeyword="annot",fileBase="",outDir="/Users/adam/test/", annotDf, annotColName,foundColword,verbose=TRUE){
  #some local helper functions to print debug reports, create outfile names
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner = function(x)paste(paste("lncRNA subsection= ",foundColword),x,sep="\n")
  
  
  # prime the df's gene_id
  #lncDf$gene_id_short <- apply(lncDf,1,function(x)as.vector(strsplit(x[c("gene_id")],"\\."))[[1]][1])
  
  # make which genes have been found in annotation
  lncDf$withinSubset = "false"
  lncDf[which(lncDf$gene_id_short %in% annotDf$gene_id_short),]$withinSubset = "true"
  
  ddply(annotDf,.(gene_id),function(x)x[1]) -> ensLnc.df
  getLncName <- function(x){y<-ensLnc.df[which(ensLnc.df$gene_id == x ), "lncRnaName"][1];if(is.na(y)){"notFound"}else{y}}
  lncDf[[annotColName]] <- as.vector(sapply(lncDf[["gene_id_short"]],function(x)getLncName(x)))
  #lncDf[is.na(lncDf[[annotColName]]),annotColName] <- "notFound"
  
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
  
  #p#am.pam.melt <- melt(pap.pam.df[c(".id","exprInBoth", "none","polyAOnly","nonPolyAOnly")],id.vars=".id")
   #p#am.pam.melt$variable <- factor(pam.pam.melt$variable,levels=c("exprInBoth","polyAOnly","nonPolyAOnly", "none"))
   #lnpa.df<-df[c("gene_id",expr.cols.nonPolyA)]
   #lpa.df<-df[c("gene_id",expr.cols.polyA)]
   #colnames(lpa.df) <- c("gene_id",expr.cells)
   #colnames(lnpa.df) <- c("gene_id",expr.cells)
#   
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
  
  
  # okay, now some plots: 
  ###### BEGIN PCA
#   lnc.pca         <- princomp(lncDf[exprCols])
#   lnc.pca.factors <-  as.data.frame(lnc.pca$scores)
#   lnc.pca.factors[["gene_id"]] <- lncDf[["gene_id"]]
#   lnc.pca.factors[["withinSubset"]] <-lncDf[["withinSubset"]]
#   lnc.pca.factors[["lncRnaName"]] <-lncDf[["lncRnaName"]]
#   
#   lnc.pca.factors[["Comp.1"]] <-lnc.pca.factors[["Comp.1"]] - min(lnc.pca.factors[["Comp.1"]]) + 1
#   lnc.pca.factors[["Comp.2"]] <-lnc.pca.factors[["Comp.2"]] - min(lnc.pca.factors[["Comp.2"]]) + 1
#   PCA.title <- titleWithBanner( "PCA of [expression vectors] foreach lnc ")
#  ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
#    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
#     layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
#      ggtitle(PCA.title) +
#       xlab("log(First Component)") + ylab("log(Second Component)")
#   ggsave(file=makeOutFile("lncRNA-pca-exprVectors.pdf"),height=8,width=8)
#   
#   lnc.sub <- which(lnc.pca.factors$withinSubset == "true")
#   lnc.min <- lnc.pca.factors[which(lnc.pca.factors[["Comp.1"]] != min(lnc.pca.factors[["Comp.1"]])),]
#   
#   ggplot(lnc.min,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
#     layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
#     layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
#    ggtitle(PCA.title)+
#     xlab("log(First Component)")+
#    ylab("log(Second Component)")+
#    ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min.pdf"),height=8,width=8)
  
  ####### END PCA
  
  
  
  
  
  
  
  ##############################  ##############################  ##############################  ##############################
  ##############################  ##############################  ##############################  ##############################
  ##############################  ##############################  ##############################  ##############################
  ##############################  ##############################  ##############################  ##############################
  # Experimental plots...

  
  
  ##############################  ##############################  ##############################  ##############################  ##############################  ##############################  ##############################  ##############################
  ##############################  ##############################  ##############################  ##############################
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
  
  
  
  
dispersion.df <- merge(lncDf,combinedDerr.df[c("gene_id",der.cols)],by="gene_id")
  
  
  disp.cols <- c(der.cols,measure)
  
  
makeDispersionPlot <- function(plot.df,colName){
 outSapplyLocal <- paste("lncRNA-",colName,"-dispersion-plot-subSetLabelled.pdf",sep="")
 title <- paste("jitter plot of ",colName," data\n","subsetted lncs w/ label",sep="")
 x.high = max(log(plot.df[[colName]])+2)
 x.low =  min(log(plot.df[[colName]])+4)
 
 plot.df[["hjust"]] <- runif(dim(plot.df)[1]) -0.5
 #    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
 #    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
 colFunStr <- paste("log(",colName,")",sep="")
 found.df <- plot.df[plot.df$lncRnaName != "notFound",]
 found.df[which(is.finite(log(found.df[[colName]])) != 0),]
 
 notFound.df <- plot.df[plot.df$lncRnaName == "notFound",]
 print("intial munge")
 labels.df <- getHjustForRange(log(as.numeric(found.df[[colName]])),found.df[["lncRnaName"]])
 str(labels.df)
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

} #

calcOffset <-function(x){if(runif(1) < 0.5){ seq(0,1000,by=0.12)[seq_along(x$vValue)]}else{ 1-seq(0,1000,by=0.12)[seq_along(x$vValue)]}}

getHjustForRange <- function(values,names){
  notInf <- which(is.finite(values)) 
  values <- values[notInf]
  names <- names[notInf]
  
  col.min <-min(values)
   col.max <- max(values)
   print(c(col.min,col.max))
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
  
  
}



pcaAnalysis <- function(lncDf,annotDf){
  
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner = function(x)paste(paste("lncRNA subsection= ",foundColword),x,sep="\n")
  
  
  # prime the df's gene_id
  #lncDf$gene_id_short <- apply(lncDf,1,function(x)as.vector(strsplit(x[c("gene_id")],"\\."))[[1]][1])
  
  # make which genes have been found in annotation
  lncDf$withinSubset = "false"
  lncDf[which(lncDf$gene_id_short %in% annotDf$gene_id_short),]$withinSubset = "true"
  
  ddply(annotDf,.(gene_id),function(x)x[1]) -> ensLnc.df
  getLncName <- function(x){ensLnc.df[which(ensLnc.df$gene_id == x ), "lncRnaName"][1]}
  lncDf[c(annotColName)] <- apply(lncDf,1,function(x)getLncName(x[c("gene_id_short")]))
  lncDf[is.na(lncDf[c(annotColName)]),annotColName] <- "notFound"
  
  lnc.pca         <- princomp(lncDf[exprCols])
  lnc.pca.factors <-  as.data.frame(lnc.pca$scores)
  lnc.pca.factors[["gene_id"]] <- lncDf[["gene_id"]]
  lnc.pca.factors[["withinSubset"]] <-lncDf[["withinSubset"]]
  lnc.pca.factors[["lncRnaName"]] <-lncDf[["lncRnaName"]]
  
  lnc.pca.factors[["Comp.1"]] <-lnc.pca.factors[["Comp.1"]] - min(lnc.pca.factors[["Comp.1"]]) + 1
  lnc.pca.factors[["Comp.2"]] <-lnc.pca.factors[["Comp.2"]] - min(lnc.pca.factors[["Comp.2"]]) + 1
  PCA.title <- titleWithBanner( "PCA of [expression vectors] foreach lnc ")
  ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(PCA.title) +
    xlab("log(First Component)") + ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors.pdf"),height=8,width=8)
  
  lnc.sub <- which(lnc.pca.factors$withinSubset == "true")
  lnc.min <- lnc.pca.factors[which(lnc.pca.factors[["Comp.1"]] != min(lnc.pca.factors[["Comp.1"]])),]
  
  ggplot(lnc.min,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(PCA.title)+
    xlab("log(First Component)")+
    ylab("log(Second Component)")+
    ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min.pdf"),height=8,width=8)
}








