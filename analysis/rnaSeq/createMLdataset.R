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
  

  func.df <- createEnsList()
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











pcaAnalysis <- function(lncDf,exprCols,normalFun=FALSE,outDir,fileBase,foundColword){
  
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner = function(x)paste(paste("lncRNA subsection= ",foundColword),x,sep="\n")
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
  
  #create a heatmap of loadings
  as.data.frame(as.matrix(lnc.pca$loadings)[])-> loadings.df
  loadings.df[["cellType"]] = rownames(loadings.df)
  loadings.melt =  melt(loadings.df)

  ggplot(comp.var.melt,aes(x=component,y=value))+geom_line()+theme_bw()+facet_wrap(~variable)+
    ggtitle("PCA lncRNA expressionmatrix\nComponent variation over all components")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-egeinvals.pdf"),height=8,width=8)
  
  ggplot(data.frame(cumVar=cumsum(comp.var$pca.std.dev^2)/sum(comp.var$pca.std.dev^2),component=1:32),aes(x=component,y=cumVar))+geom_line()+geom_point()+theme_bw()+
    ggtitle("PCA lncRNA expressionmatrix\nComponent variation over all components")+xlab("cummulative component")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-egeinvals-ecdf.pdf"),height=8,width=8)
  
  ggplot(data.frame(cumVar=cumsum(comp.var$cor.std.dev^2)/sum(comp.var$cor.std.dev^2),component=1:32),aes(x=component,y=cumVar))+geom_line()+geom_point()+theme_bw()+
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
  
#   lnc.pca.factors.cov <-  as.data.frame(lnc.pca.cov$scores)
#   lnc.pca.factors.cov[["gene_id"]] <- lncDf[["gene_id"]]
#   lnc.pca.factors.cov[["withinSubset"]] <-lncDf[["withinSubset"]]
#   lnc.pca.factors.cov[["lncRnaName"]] <-lncDf[["lncRnaName"]]
#   
#   lnc.pca.factors.cov[["Comp.1"]] <-lnc.pca.factors.cov[["Comp.1"]] - min(lnc.pca.factors.cov[["Comp.1"]]) + 1
#   lnc.pca.factors.cov[["Comp.2"]] <-lnc.pca.factors.cov[["Comp.2"]] - min(lnc.pca.factors.cov[["Comp.2"]]) + 1
#   lnc.pca.factors.cov[["Comp.3"]] <-lnc.pca.factors.cov[["Comp.3"]] - min(lnc.pca.factors.cov[["Comp.3"]]) + 1
#   lnc.pca.factors.cov[["Comp.4"]] <-lnc.pca.factors.cov[["Comp.4"]] - min(lnc.pca.factors.cov[["Comp.4"]]) + 1
#   
  
  PCA.title <- titleWithBanner( "PCA of [expression vectors] foreach lnc ")
  ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(PCA.title) +
    xlab("log(First Component)") + ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors.pdf"),height=8,width=8)
  
#   PCA.title <- titleWithBanner( "PCA of COVARIANCE [expression vectors] foreach lnc ")
#   ggplot(lnc.pca.factors.cov,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
#     layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
#     layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
#     ggtitle(PCA.title) +
#     xlab("log(First Component)") + ylab("log(Second Component)")
#   ggsave(file=makeOutFile("lncRNA-pca-exprVectors-cov.pdf"),height=8,width=8)
#   
  PCA.title.cor <- titleWithBanner( "PCA of CORRELATION [expression vectors] foreach lnc ")
  ggplot(lnc.pca.factors.cor,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(PCA.title.cor) +
    xlab("log(First Component)") + ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-cor.pdf"),height=8,width=8)

  lnc.sub <- which(lnc.pca.factors$withinSubset == "true")
  lnc.min <- lnc.pca.factors[which(lnc.pca.factors[["Comp.1"]] != min(lnc.pca.factors[["Comp.1"]])),]
  
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
    ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min-comp1-3.pdf"),height=8,width=8)
  
  ggplot(lnc.min,aes(x=log(Comp.2),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title,"\nwithout outlier"))+
    xlab("log(First Component)")+
    ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-exprVectors-remove-min-comp2-3.pdf"),height=8,width=8)
  
  lnc.sub.cor <- which(lnc.pca.factors.cor$withinSubset == "true")
  lnc.min.cor <- lnc.pca.factors.cor[which(lnc.pca.factors.cor[["Comp.1"]] != min(lnc.pca.factors.cor[["Comp.1"]])),]
  
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
    ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-cor-exprVectors-remove-min-com1-3.pdf"),height=8,width=8)
  
  ggplot(lnc.min.cor,aes(x=log(Comp.2),y=log(Comp.3),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors.cor[lnc.pca.factors.cor$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)+
    ggtitle(paste(PCA.title.cor,"\nwithout outlier CORRELATION"))+
    xlab("log(First Component)")+
    ylab("log(Second Component)")
  ggsave(file=makeOutFile("lncRNA-pca-cor-exprVectors-remove-min-comp2-3.pdf"),height=8,width=8)  
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
  pcaAnalysis(lncDf=lncList.df,exprCols=lnc.colIndex,normalFun=FALSE,foundColword="foundInLiterature",outDir="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaDbPlusLit",fileBase="PCA-")
}

getLncRnaName <- function(label,ens,df){ifdf[ens,"external_gene_id"]}

getLncRNAExprDf <- function(){
  outfile <- getFullPath("data/lncRnaExpr_ml.tab")
  if (!file.exists(outfile)){
  # get all the functional or annotated lncRNAs
  r.df <- createEnsList()
  # get derrien's lnc'd list
  combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
  lnc.out.dir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/topTransForGenePlots_eachSample"
  lnc.fileBase = "lncRNA-gene=MaxTransSample"
  lnc.transKeyword = "lncRNA genes as max trans in sample"
  df.1=readInTable(lnc.in.file)
  df.1$gene_id_short <- sapply(df.1$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
  df.1[["label"]] <- 0
  df.1[which(df.1$gene_id_short %in% r.df$ensembl_gene_id),"label"] <- 1
  
  label.df <- df.1[df.1$label == 1,]
  label.df <- unique(merge(label.df,unique(r.df[c("ensembl_gene_id","external_gene_id")]),by.x="gene_id_short",by.y="ensembl_gene_id"))
  label.df$lncRnaName <- label.df$external_gene_id 
  label.df$external_gene_id <- NULL
                     
  nolabel.df <- df.1[df.1$label == 0,]
  nolabel.df$lncRnaName = nolabel.df$gene_id_short
                     
  combined.df <- rbind(label.df,nolabel.df)
  combined.expr.index = 3:34
                     
  derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  derrien.mean.df <- ddply(derrien.df,.(LncRNA_GeneId),numcolwise(mean))
                     
                     
  within(derrien.mean.df,{gene_id_short = sapply(LncRNA_GeneId,function(x)unlist(strsplit(x,"\\."))[[1]])})  -> derrien.mean.df
  combined.df <- merge(combined.df,derrien.mean.df[c("gene_id_short","LncRNA_tx_size")],by.x="gene_id_short",by.y="gene_id_short")
  exportAsTable(df=combined.df,file=outfile)                
  }                  
   else{
     readInTable(outfile)
     
   }                  
}

runEigenRankOnExprData <- function(distfile=getFullPath("./data/exprDistance.tab")){
 df <- readInTable(getFullPath("data/lncRnaExpr_ml.tab"))
 expr.mat <- as.matrix(df[3:34])
 distance.mat <- as.matrix(dist(expr.mat)) 
 d = calcEigenRank(mat=distance.mat,label=df$label[1:500])
 ggplot( d,aes(x=log(abs(rank))^(1),color=factor(type),fill=factor(type)))+geom_density(alpha=I(0.6))+theme_bw() +
   ggtitle("EigenRank of different lncRNA subset\nClusterd on cshl rna-seq expression data")
 ggsave(getFullPath("plots/eigenRank-allExprData.pdf"),width=14,height=14)
}

runEigenRankOnExprData <- function(distfile=getFullPath("./data/exprDistance.tab")){
  df <- readInTable(getFullPath("data/lncRnaExpr_ml.tab"))
  feature.mat <- as.matrix(df[1:500,c(colnames(df)[3:34],"tissSpec","LncRNA_tx_size")])
 
  distance.mat <- as.matrix(dist(feature.mat)) 
  d = calcEigenRank(mat=distance.mat,label=df$label)
  ggplot( d,aes(x=log(abs(rank))^(1),color=factor(type),fill=factor(type)))+geom_density(alpha=I(0.6))+theme_bw() + 
    ggtitle("EigenRank of lncRNA subsets \nClustered on transcript length & tissue specifity + expression data")
  ggsave(getFullPath("plots/eigenRank-allExprData+tissSpec+translength.pdf"),width=14,height=14)
}



calcEigenRankDf <- function(lab,nolab,cols=3:34,scaleFun=function(x)x){
    
  #mat.df = scaleFun(rbind(lab[,cols],nolab[,cols]))
  mat.df = rbind(lab[,cols],nolab[,cols])
  
  mat = as.matrix(dist(mat.df))
  mat = mat^(-1)
  diag(mat) = 0
  mat[which(mat == Inf)] = 0
  #mat %*% eigen(x=mat)$vector[,1] -> x1
  df <- data.frame(rank=(eigen(x=mat)$vector[,1]),label = c(lab$label,nolab$label))
  #data.frame(label = c(lab$label,nolab$label), lncRnaName = c(lab$lncRnaName,nolab$lncRnaName))
  df$rank = (df$rank / sum(df$rank)) * length(df$rank)
  df
}
calcEigenRankLabelNoLabel <- function(lab,nolab,cols=3:34,scaleFun=function(x)x){
  
  #mat.df = scaleFun(rbind(lab[,cols],nolab[,cols]))
  mat.df = rbind(lab[,cols],nolab[,cols])
  mat = as.matrix(dist(mat.df))
  df <- calcEigenRankIter(mat=mat,label = c(lab$label,nolab$label))
  df$lncRnaName = c(lab$lncRnaName,nolab$lncRnaName)
  df
}

calcEigenRankLabelNoLabelSwitch <- function(lab,nolab,cols=3:34,scaleFun=function(x)x){
  
  #mat.df = scaleFun(rbind(lab[,cols],nolab[,cols]))
  mat.df = rbind(nolab[,cols],lab[,cols])
  matFromDf = as.matrix(mat.df)
  mat = as.matrix(dist(mat.df))
  
  df <- calcEigenRankIter2(mat=mat,label = c(nolab$label,lab$label),norm=FALSE)
  df$lncRnaName = c(nolab$lncRnaName,lab$lncRnaName)
  df
}


plotAndCalcEigenRank <- function(df,
                                 outFile=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr.pdf"),
                                 title="EigenRank of lncRNA subsets \nfeatures are expression data\nlabel: 1=functional, 0=unlabled ",
                                cols=2:33)
{
label.df = df[which(df$label == 1),] 
nolabel.df = df[which(df$label == 0),] 

nolabel.df$index = seq_along(nolabel.df$label)
nolabel.df$group = cut(nolabel.df$index,pretty(nolabel.df$index,10))

eigen.df=ddply(nolabel.df,.(group),function(x) calcEigenRankDf(lab = label.df, nolab = x,cols=cols))
#eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval = wilcox.test(x[x$label == 0,"eigenVecValue"],x[ x$label == 1,"eigenVecValue"])))

ggplot(eigen.df, aes( x = log(rank), fill = factor(label))) + geom_density(alpha=I(0.6))+facet_wrap( ~group) + theme_bw() + 
  ggtitle(title)
ggsave(file=outFile)

}





eigenRankSplitUp <- function(folder = getFullPath("plots/rnaSeq-eigenRank")){
  
  df <- readInTable(getFullPath("data/lncRnaExpr_ml.tab"))
  lnpa <- grep("longNonPolyA$",colnames(df))
  lpa <- grep("longPolyA$",colnames(df))
  
  double.cols = as.vector(which(unlist(sapply(colnames(df), function(x)(typeof(df[[x]]) == "double")))))
  df.norm <- as.data.frame(sapply(colnames(df[,double.cols]),function(x){(df[[x]] - mean(df[[x]])) /sd(df[[x]])}))
  df.norm[["label"]] = df[["label"]]
  df.norm[["lncRnaName"]] = df[["lncRnaName"]]
  df.norm.cols = as.vector(which(unlist(sapply(colnames(df.norm), function(x)(typeof(df.norm[[x]]) == "double")))))
  
  plotAndCalcEigenRank(df=df.norm,
                       outFile=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr_normAllCols.pdf"),
                                                       title="EigenRank of lncRNA subsets \nfeatures are expression data + statistics(normalized)\nlabel: 1=functional, 0=unlabled ",
                                                       cols=df.norm.cols)
  
  plotAndCalcEigenRank(df=df.norm,
                       outFile=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr_longNonPolyA.pdf"),
                       title="EigenRank of lncRNA subsets \nLongNonPolyA RNA seq(normalized)\nlabel: 1=functional, 0=unlabled ",
                       cols= grep("longNonPolyA$",colnames(df)))
  
  plotAndCalcEigenRank(df=df.norm,
                       outFile=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr_longPolyA.pdf"),
                       title="EigenRank of lncRNA subsets \nlongPolyA RNA seq(normalized)\nlabel: 1=functional, 0=unlabled ",
                       cols= grep("longPolyA$",colnames(df.norm)))
  
  plotAndCalcEigenRank(df=df.norm,
                       outFile=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr_exprNorm.pdf"),
                       title="EigenRank of lncRNA subsets \nAll RNA seq(lpa + lnpa)(normalized)\nlabel: 1=functional, 0=unlabled ",
                       cols= c(grep("longNonPolyA$",colnames(df)),grep("longPolyA$",colnames(df.norm))))
  

  
  
  
  
  
  label.df = df[which(df$label == 1),] 
  nolabel.df = df[which(df$label == 0),] 
  nolabel.df$index = seq_along(nolabel.df$label)
  nolabel.df$group = cut(nolabel.df$index,pretty(nolabel.df$index,10))
  eigen.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x))
  #eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval=wilcox.test(x[x$label == 0,"eigenVecValue"],x[x$label == 1,"eigenVecValue"])))
  ggplot(eigen.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nfeatures are expression data\nlabel: 1=functional, 0=unlabled ")

  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr.pdf"))
  
  
  eigen1.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x,cols=c("LncRNA_tx_size","tissSpec","averageExpr","entropyExpr")))
  #eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval=wilcox.test(x[x$label == 0,"eigenVecValue"],x[x$label == 1,"eigenVecValue"])))
  ggplot(eigen1.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nfeatures are transcript size, aveExp,tissSpec,entropy")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubset-txSize-tissSpec-aveExpr-entropy.pdf"))
  
  eigen1.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x,cols=c(colnames(), "LncRNA_tx_size","tissSpec","averageExpr","entropyExpr")))
  #eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval=wilcox.test(x[x$label == 0,"eigenVecValue"],x[x$label == 1,"eigenVecValue"])))
  ggplot(eigen1.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nfeatures are transcript size, aveExp,tissSpec,entropy")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubset-txSize-tissSpec-aveExpr-entropy.pdf"))
  
  #try with 
  eigen.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x,cols=c("LncRNA_tx_size","averageExpr"),scaleFun=numcolwise(mean)))
  ggplot(eigen.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nClustered on transcript length & tissue specifity + expression data")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubset-aveExpr-transLen.pdf"))
}

eigenRankSplitUpIteration <- function(folder = getFullPath("plots/rnaSeq-eigenRank")){
  df <- readInTable(getFullPath("data/lncRnaExpr_ml.tab"))
  label.df=df[which(df$label == 1),] 
  nolabel.df=df[which(df$label == 0),] 
  nolabel.df$index = seq_along(nolabel.df$label)
  nolabel.df$group = cut(nolabel.df$index,pretty(nolabel.df$index,10))
  eigen.df=ddply(nolabel.df,.(group),function(x)calcEigenRankLabelNoLabel(lab=label.df,nolab=x))
  #eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval=wilcox.test(x[x$label == 0,"eigenVecValue"],x[x$label == 1,"eigenVecValue"])))
  ggplot(eigen.df,aes(x=log(rank),fill=factor(type)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nfeatures are expression data\nlabel: 1=functional, 0=unlabled ")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubsetExpr-iteration.pdf"))
  
  eigen1.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x,cols=c("LncRNA_tx_size","tissSpec","averageExpr","entropyExpr")))
  #eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval=wilcox.test(x[x$label == 0,"eigenVecValue"],x[x$label == 1,"eigenVecValue"])))
  ggplot(eigen1.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nfeatures are transcript size, aveExp,tissSpec,entropy")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubset-txSize-tissSpec-aveExpr-entropy.pdf"))
  
  eigen1.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x,cols=c(colnames(), "LncRNA_tx_size","tissSpec","averageExpr","entropyExpr")))
  #eigenpvals.df=ddply(nolabel.df,.(group),function(x)data.frame(pval=wilcox.test(x[x$label == 0,"eigenVecValue"],x[x$label == 1,"eigenVecValue"])))
  ggplot(eigen1.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nfeatures are transcript size, aveExp,tissSpec,entropy")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubset-txSize-tissSpec-aveExpr-entropy-iter.pdf"))
  
  #try with 
  eigen.df=ddply(nolabel.df,.(group),function(x)calcEigenRankDf(lab=label.df,nolab=x,cols=c("LncRNA_tx_size","averageExpr"),scaleFun=numcolwise(mean)))
  ggplot(eigen.df,aes(x=log(abs(eigenVecValue))^(1),fill=factor(label)))+geom_density(alpha=I(0.6))+facet_wrap(~group)+theme_bw()+ 
    ggtitle("EigenRank of lncRNA subsets \nClustered on transcript length & tissue specifity + expression data")
  ggsave(file=getFullPath("plots/rnaSeq-eigenRank/eigenRanksBySubset-aveExpr-transLen.pdf"))
}

createDatasetOfPcaPlusLength <- function(df){
  
  combined.df <- getLncRNAExprDf()
  
  lnc.pca             <- princomp(combined.df[combined.expr.index],cor=FALSE)
  lnc.pca.factors <-  as.data.frame(lnc.pca$scores)
  lnc.pca.factors <- cbind(combined.df[c("gene_id_short","label","lncRnaName","LncRNA_tx_size")],lnc.pca.factors)
 
 ## lnc.pca.factors[["Comp.1"]] <-lnc.pca.factors[["Comp.1"]] - min(lnc.pca.factors[["Comp.1"]]) + 1
  #lnc.pca.factors[["Comp.2"]] <-lnc.pca.factors[["Comp.2"]] - min(lnc.pca.factors[["Comp.2"]]) + 1
  ##lnc.pca.factors[["Comp.3"]] <-lnc.pca.factors[["Comp.3"]] - min(lnc.pca.factors[["Comp.3"]]) + 1
  #lnc.pca.factors[["Comp.4"]] <-lnc.pca.factors[["Comp.4"]] - min(lnc.pca.factors[["Comp.4"]]) + 1
  
   exportAsTable(df=lnc.pca.factors[c("gene_id_short","lncRnaName","label","LncRNA_tx_size","Comp.1","Comp.2","Comp.3")],file=outfile)
                     lnc.pca.factors[c("gene_id_short","lncRnaName","label","LncRNA_tx_size","Comp.1","Comp.2","Comp.3")]                     
}


getMLdataset <- function(){
  infile <- getFullPath("data/lncRnaPcaTransLength.tab")
  if (!file.exists(infile)){
    createDatasetOfPcaPlusLength()
  }
  else {
    readInTable(infile)
  }
}


calcEigenRank <- function(mat=as.matrix(dist(iris[,2:4])),label=iris$Species){
mat = mat^(-1)
diag(mat) = 0
mat[which(mat == Inf)] = 0
mat %*% eigen(x=mat)$vector[,1] -> x1

  vec <- x1;
#for (i in 1:100){
#  vec.tmp = mat %*% (vec/sum(vec))
##  vec = vec.tmp
#}
  
  df <- data.frame(rank=vec,type=label)
  df
}

calcEigenRankIter <- function(mat=as.matrix(dist(iris[,2:4])),label=iris$Species,norm=TRUE,diag=0,printIteration=FALSE){
  mat = mat^(-1)
  diag(mat) = diag
  mat[which(mat == Inf)] = 0
  mat %*% eigen(x=mat)$vector[,1] -> x1
  #vec = runif(length(label))
  vec <- x1;
  i = 0;
 repeat {
   i = i + 1;
    vec.tmp = mat %*% (vec/sum(vec))
   
   if(printIteration == TRUE){ 
   print(paste("iteration",i,"threshold = ",sum(vec - vec.tmp)))
   } 
   if(abs(sum(vec - vec.tmp)) < 10^-10 || i > 4000){
      break;
    }
    vec = vec.tmp
  }
 
  
  df <- data.frame(rank=vec,type=label)
  if (norm == TRUE){
     df$rank = (df$rank / sum(df$rank)) * length(df$rank)  ## (value / sum) * count = average value is one
     df
  }else{
    df
  }
  
}
calcEigenRankIter2 <- function(mat=as.matrix(dist(iris[,2:4])),label=iris$Species,norm=TRUE,diag=0,printIteration=FALSE){
  mat = mat^(-1)
  diag(mat) = diag
  mat[which(mat == Inf)] = 0
  mat %*% eigen(x=mat)$vector[,1] -> x1
  
  vec <- x1;
  i = 0;
 
  matR <<- mat
  #mat = matC
  print("got here")
  repeat {
    i = i + 1;
    vec.tmp = mat %*% (vec/sum(vec))
    
    if(printIteration == TRUE){ 
      print(paste("iteration",i,"threshold = ",sum(vec - vec.tmp)))
    } 
    if(abs(sum(vec - vec.tmp)) < 10^-10 || i > 4000){
      vec = vec.tmp
      break;
    }
    vec = vec.tmp
  }
  
  
  df <- data.frame(rank=vec,type=label)
  if (norm == TRUE){
    df$rank = (df$rank / sum(df$rank)) * length(df$rank)  ## (value / sum) * count = average value is one
    df
  }else{
    df
  }
  
}


testEigenRank <- function(){
  df.1 <- calcEigenRank()
  mat.df = as.matrix(dist(matrix(data=runif(3*150),ncol=3,nrow=150)))
  df.2 <- calcEigenRank(mat=mat.df,label="random")
  df.3 <- calcEigenRankIter()
  df.4 <- calcEigenRankIter(mat=mat.df,label="random")
  #create a cluster within a distribution 
  
  
  
  df.1[["facet"]] = "actualData"  
  df.2[["facet"]] = "randomUniformData"  
  
  df.3[["facet"]] = "actualData_Iter"  
  df.4[["facet"]] = "randomUniformData_Iter" 
  
  
  
  comb.df <- rbind(df.1,df.2,df.3,df.4)
  
  ggplot(comb.df,aes(x=log(rank),fill=factor(type)))+geom_bar(binwidth=0.1)+facet_wrap(~facet)
  
  
  matLabel <- c(rep(1,20),rep(0,200))
 
  clust.mat =  rbind(cbind(rnorm(20,mean=0,sd=1),rnorm(20,mean=0,sd=1)),cbind(rnorm(200,mean=0,sd=30),rnorm(200,mean=0,sd=30)))
  matClust <- as.matrix(dist(clust.mat))
  clust.df= as.data.frame(clust.mat)
  colnames(clust.df) <- c("x","y")
  clust.df[["label"]] = matLabel
  clust.df[["type"]] = "clustered"
  
  rand.mat = rbind(cbind(rnorm(20,mean=0,sd=30),rnorm(20,mean=0,sd=30)),cbind(rnorm(200,mean=0,sd=30),rnorm(200,mean=0,sd=30)))
  matRand <-  as.matrix(dist(rand.mat))
  rand.df = as.matrix(dist(rand.mat))
  rand.df = as.data.frame(rand.mat)
  colnames(rand.df) = c("x","y")
  rand.df[["label"]] = matLabel
  rand.df[["type"]] = "noCluster"
  
  
  df.clust.iter <- calcEigenRankIter(mat=matClust,label=matLabel,diag=0)
  df.clust <- calcEigenRank(mat=matClust,label=matLabel)
  df.rand.iter <- calcEigenRankIter(mat=matRand,label=matLabel,diag=0)
  df.rand <- calcEigenRank(mat=matRand,label=matLabel)
  df.clust[["facet"]] = "clusteredData"  
  df.rand[["facet"]] = "randomData"  
  df.clust.iter[["facet"]] = "clusteredDataIter"  
  df.rand.iter[["facet"]] = "randomDataIter"  
  
  comb.df <- rbind(df.clust,df.rand,df.clust.iter,df.rand.iter)
  comb.df= ddply(comb.df,.(facet),transform,norm.EigenVector=sum(rank))
  comb.df= ddply(comb.df,.(facet),transform,count.EigenVector=length(rank))
  
  
  ggplot(rbind(df.rand.iter,df.clust.iter),aes(x=as.numeric(rank),fill=factor(type)))+geom_bar()+facet_wrap(~facet)+theme_bw()+scale_fill_manual(values=c("pink","black"))+ ggtitle("pageRank on 2d points\nRight Plot: sd = 30 & sd = 1,Left Plot: sd = 30 both groups\n20 black & 200 pink ")
  ggsave("~/Desktop/matrix-eigen.pdf",height=4)
  
  ggplot(rbind(rand.df,clust.df),aes(x=x,y=y,color=factor(label),size=2))+geom_point(alpha=I(0.9),size=3)+theme_bw()+facet_wrap(~type)+scale_color_manual(values=c("pink","black"))
  ggsave("~/Desktop/matrix-data.pdf",height=4)
  
  clust.df <- data.frame(x=df.rand$rank*length(df.rand$rank)/sum(df.rand$rank),y=df.rand.iter$rank*length(df.rand.iter$rank)/sum(df.rand.iter$rank))
  ggplot(data.frame(x=df.rand$rank,y=df.rand.iter$rank),aes(x,y))+geom_point()  
  ggsave("~/Desktop/eigenVecVsIter.pdf")  
}
