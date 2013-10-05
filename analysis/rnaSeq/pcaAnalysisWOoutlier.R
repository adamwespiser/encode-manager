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


testFun <- function(){
  iris.norm <- as.data.frame(sapply(colnames(iris[,1:4]),function(x){(iris[[x]] - mean(iris[[x]])) /sd(iris[[x]])}))
  x0 = princomp(as.matrix(iris[,1:4]),cor=FALSE,scale=FALSE)
  x0$scores[,1]
  
  x01 = prcomp(as.matrix(iris[,1:4]),cor=FALSE,scale=FALSE)
  x01$rotation[,1]  # x0$loadings[,1]
  drop(scale(as.matrix(iris[,1:4]), center = x01$center, scale = x01$scale) %*% x01$rotation[,1])  # x0$scores[,1]
  ex
  x1 = princomp(as.matrix(iris.norm[,1:4]),cor=FALSE,scale=FALSE)
  x2 = princomp(as.matrix(iris[,1:4]),cor=FALSE,scale=TRUE)
  
  data("heptathlon", package = "HSAUR")
  score <- which(colnames(heptathlon) == "score")
  heptathlon_pca <- prcomp(heptathlon[, -score], scale = TRUE)
  a1 <- heptathlon_pca$rotation[,1]
  center <- heptathlon_pca$center
  scale <- heptathlon_pca$scale
  hm <- as.matrix(heptathlon[,-score])
  drop(scale(hm, center = center, scale = scale) %*% heptathlon_pca$rotation[,1])
}


testLnc <- function(){
  princomp(comb.out[3:34],cor=FALSE)-> x10
  firstFive <- (scale(head(comb.out[,3:34]),center=x10$center,scale=x10$scale) %*% x10$loadings )[,1:5]
}

pdTest <- function(){
  df <- iris
  df.center= suppressMessages(melt(ddply(df,.(),numcolwise(mean)))[["value"]])
  df[["dist"]] = apply(df[exprCols],1,function(x) sqrt(sum((x -df.center)^2)))
  df[["withinSubset"]] = ifelse(runif(150) < 0.1,1,0)
  
  
  dfo = df[order(max(df$dist) - df$dist),] # order from highest to lowest
  dfo[["cs"]] = cumsum(dfo$withinSubset)
  dfo[["predict"]] = seq_along(dfo$withinSubset)
  dfo[["totalYone"]] = sum(dfo$withinSubset)
  dfo[["totalMembers"]] = length(dfo$withinSubset)
  
  dfo.stats =within(dfo,{
  TP = cs
  FP = predict - cs
  TN = (totalMembers - predict) - (totalYone - cs)
  FN = totalYone - cs
  sens = TP / (TP + FN)
  prec = TP / (TP + FP)
  })
  
  
  ggplot(dfo.stats,aes(y=prec,x=sens,color=predict))+geom_point(size=1)+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
}





plotDistanceAwayRatio <- function(df,exprCols,
                                  df.center= suppressMessages(melt(ddply(df[exprCols],.(),numcolwise(mean)))[["value"]]),rand=FALSE){
  
  df[["withinSubset"]] = ifelse(df[["withinSubset"]] == "true",1,0)
  if(rand == TRUE){
    df[["rand"]] = runif(length(df[["withinSubset"]]))
    df[["withinSubset"]] = df[order(df[["rand"]]),"withinSubset"]
  }
  
  df[["dist"]] = apply(df[exprCols],1,function(x) sqrt(sum((x - df.center)^2)))
  dfo = df[order(max(df$dist) - df$dist),] # order from highest to lowest
  dfo[["TP"]] = cumsum(dfo$withinSubset)
  dfo[["predict"]] = seq_along(dfo$withinSubset)
  dfo[["P"]] = sum(dfo$withinSubset)
  dfo[["totalMembers"]] = length(dfo$withinSubset)
  dfo[["N"]] = dfo[["totalMembers"]] - dfo[["P"]]
  dfo.stats =within(dfo,{
    FP = predict - TP
    TN = (totalMembers - predict) - (P - TP)
    FN = P - TP
    sens = TP / (TP + FN)
    prec = TP / (TP + FP)
    TPR = TP / P
    FPR = FP / N
    FNR = FN / P
    TNR = TN / N
    Fscore = ( 2 * prec * sens ) / ( prec + sens )
    accuracy = (TP + TN) / (P + N)
    errorRate = (FP + FN) / (P + N)
  })
  dfo.stats[c("predict","totalMembers","P","N","TP","FP","TN","FN","sens","prec","TPR","FPR","FNR","TNR","Fscore","accuracy","errorRate")]
}

saveToFilePlotDistAway <- function(df,exprCols,outfile,titleMsg){
  df = plotDistanceAwayRatio(df,exprCols)
  
  ggplot( df, aes(x=sens,y=prec)) + geom_point() +theme_bw()+
    ggtitle(titleMsg)
  ggsave(paste(outfile,"circClassifierPCA.pdf", sep=""),)
  
}


# create a heat map of the loadings for each of the principal compenents
plotLoadingsToFile <- function(pca,outfile,titleMsg=""){
  
  #extract the loadings
  as.data.frame(as.matrix(pca$loadings)[])-> loadings.df
  loadings.df[["cellType"]] = rownames(loadings.df) # push rownames to col
  
  # melt into geom_raster accectable format
  loadings.melt =  melt(loadings.df,id.var="cellType")
  loadings.melt[["component"]] = sapply(loadings.melt$variable,function(x){as.vector(strsplit(as.character(x),"\\.")[[1]])[2]})
  
  # plot the heatmap with grey neutral, and green/red as plus minus loading 
  ggplot(loadings.melt,aes(y=cellType,x=as.numeric(component),fill=value))+geom_raster()+xlab("Comp") + 
    scale_fill_gradient2(low="darkred",mid="lightgrey",high="green") + theme_bw()+
    ggtitle(titleWithBanner(paste("PCA loadings :",titleMsg,sep="")))
  ggsave(file=outfile,height=8,width=8)
  
}

#plot the variance for each component, as well the as cummulative variance up to component N, for al components in PCA
plotScreeToFile <- function(pca,outfile,outfileECDF,titleMsg=""){
  comp.var =  data.frame(pca.std.dev=pca$sdev^2,
                         component=seq_along(colnames(pca$scores))) 
  comp.var.melt <- melt(comp.var, id.var="component")  
  
  ggplot(comp.var.melt,aes(x=component,y=value))+geom_line()+theme_bw()+
    ggtitle("PCA lncRNA expressionmatrix\nComponent variation over all components")
  ggsave(file=outfile,height=8,width=8)
  
  ggplot(data.frame(cumVar=cumsum(comp.var$pca.std.dev^2)/sum(comp.var$pca.std.dev^2),component=seq_along(colnames(pca$scores))),aes(x=component,y=cumVar))+geom_line()+geom_point()+theme_bw()+
    ggtitle(titleWithBanner("variance of PCA components"))+xlab("cummulative component")
  ggsave(file=outfileECDF,height=8,width=8)
  
  
}

#take a pca object and scale the first 4 componets s.t. they are acceptable to log scaling
normalizePcaFactors <- function(pca,lncDataFrame){
  lnc.pca.factors <-  as.data.frame(pca$scores)
  within(as.data.frame(pca$scores),{
    gene_id =  lncDataFrame[["gene_id"]]
    withinSubset = lncDataFrame[["withinSubset"]]
    lncRnaName =  lncDataFrame[["lncRnaName"]]
    Comp.1 = Comp.1 - min(Comp.1) + 1
    Comp.2 = Comp.2 - min(Comp.2) + 1
    Comp.3 = Comp.3 - min(Comp.3) + 1
    Comp.4 = Comp.4 - min(Comp.4) + 1})
}

#pca 
transformDataByPCA <- function(df,
                               df.center= suppressMessages(melt(ddply(df,.(),numcolwise(mean)))[["value"]]),
                               df.scale=rep(1,length(colnames(df))),
                               pca){
  
  scale(df,center=df.center,scale=df.scale) %*% pca$loadings
}

normalizePcaFactorsFromLoading <- function(pca,lncDataFrame,exprCols){
  lncDataFrame[exprCols] = transformDataByPCA(df=lncDataFrame[,exprCols],pca=pca)
  
  if (is.numeric(exprCols)){
  colnames(lncDataFrame)[exprCols] = paste0("Comp.",seq_along(exprCols))
  } else {
  colnames(lncDataFrame)[which(colnames(lncDataFrame) %in% exprCols)] = paste0("Comp.", seq_along(exprCols))
  }
  
  normLncDf <- apply(lncDataFrame[ paste0("Comp.", seq_along(exprCols))],2,function(x)(x - min(x))+1)
  
  df.out <- data.frame(
    gene_id =  lncDataFrame[["gene_id"]],
    withinSubset = lncDataFrame[["withinSubset"]],
    lncRnaName =  lncDataFrame[["lncRnaName"]])
  
  cbind(df.out,normLncDf)

  
}

normalizePcaFactorsFromLoadingForRatioTest <- function(pca,lncDataFrame,exprCols){
  lncDataFrame[exprCols] = transformDataByPCA(df=lncDataFrame[,exprCols],pca=pca)
  colnames(lncDataFrame)[exprCols] = paste0("Comp.",seq_along(exprCols))
  lncDataFrame
}
# plot the comparisons between the first 3 components to file
plotPcaCompToFile <- function(pcaFactors,baseDir,fileNameMsg,titleMsg,comp1range=c(-Inf,Inf),comp2range=c(-Inf,Inf),comp3range=c(-Inf,Inf),labelsSet=FALSE,ptSize=I(0.2)){
PCA.title = paste("PCA of lncRNA expression data",titleMsg,sep="\n")
  
#create filenames to reflect the component comparisons
oneTwo = paste(baseDir,"comps-1-2-",fileNameMsg,".pdf",sep="")  
oneThree = paste(baseDir,"comps-1-3-",fileNameMsg,".pdf",sep="")  
twoThree = paste(baseDir,"comps-2-3-",fileNameMsg,".pdf",sep="")  

#pcaFactors[pcaFactors$Comp.1 == min(pcaFactors$Comp.1) || pcaFactors$Comp.1 == min(pcaFactors$Comp.1), "withinSubset"] = true
#pcaFactors[pcaFactors$Comp.2 == min(pcaFactors$Comp.2) || pcaFactors$Comp.2 == min(pcaFactors$Comp.2), "withinSubset"] = true
#pcaFactors[pcaFactors$Comp.3 == min(pcaFactors$Comp.3) || pcaFactors$Comp.3 == min(pcaFactors$Comp.3), "withinSubset"] = true

pcaFactors = within(pcaFactors,{
  logComp1  = log(Comp.1)
  logComp2  = log(Comp.2)
  logComp3  = log(Comp.3)
  range1    = sapply(logComp1, function(x) withinRange3(x = x, comp1range[1], comp1range[2]))
  range2    = sapply(logComp2, function(x) withinRange3(x = x, comp2range[1], comp2range[2]))
  range3    = sapply(logComp3, function(x) withinRange3(x = x, comp3range[1], comp3range[2]))
  func      = withinSubset == "true"
  comp12set = (range1 | range2) | func
  comp13set = (range1 | range3) | func
  comp23set = (range2 | range3) | func 
  lncRnaName = as.vector(sapply(lncRnaName,function(x)sub(x=x,"ENSG00000","")))
  functionalLncRna = withinSubset
  })

#pcaFactor[pcaFactors$logComp1 == min(pcaFactors$logComp1) || pcaFactors$logComp1 == min(pcaFactors$logComp1), ]
if (1 == 0){
oneTwoA = paste(baseDir,"comps-1-2-A",fileNameMsg,".pdf",sep="")  

ggplot(pcaFactors,aes(x = logComp1, y = logComp2, color = withinSubset,fill = withinSubset))+stat_density2d(alpha=I(0.3))+
  theme_bw()+stat_density2d(data=pcaFactors,alpha=I(0.1),geom="polygon")
ggsave(file=oneTwoA,height=8,width=8)

# ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width,fill=Species))+stat_density2d(alpha=I(0.1),geom="polygon")+theme_bw()
oneTwoB = paste(baseDir,"comps-1-2-B",fileNameMsg,".pdf",sep="")  

ggplot(pcaFactors,aes(x = logComp1, y = logComp2,fill = withinSubset))+stat_density2d(n=1000,alpha=I(0.1),geom="polygon",countour=FALSE)+theme_bw()
ggsave(file=oneTwoB,height=8,width=8)
} 

#ggplot of lncRNA expression component factors  Two additional layers are added for the annotated points
g1 <- ggplot(pcaFactors, aes(x = logComp1, y = logComp2, color = withinSubset, shape = withinSubset)) + geom_point( size =ptSize) + theme_bw() + 
 # layer(data = pcaFactors[pcaFactors$comp12set ,], mapping = aes(x = logComp1, y = logComp2, hjust = 0.9,vjust=-1,label=lncRnaName),geom = "text",size = 4) +
  layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp1, y = logComp2),size =  I(ptSize), geom = "point") +
#  layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp1, y = logComp2, hjust = 0.9,vjust=-1,label=sub(x=lncRnaName,pattern="ENSG00000",replacement="")),geom = "text",size = 4) +
  ggtitle(PCA.title) + 
  xlab("log(First Component)") +
  ylab("log(Second Component)")
if (labelsSet == TRUE){
  g1 = g1 + layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp1, y = logComp2, hjust = 0.9,vjust=-1,label=sub(x=lncRnaName,pattern="ENSG00000",replacement="")),geom = "text",size = 4) +    
  layer(data = pcaFactors[pcaFactors$comp12set ,], mapping = aes(x = logComp1, y = logComp2, hjust = 0.9,vjust=-1,label=lncRnaName),geom = "text",size = 4)
}

g1;ggsave(file=oneTwo,height=8,width=8)

g2 <- ggplot(pcaFactors, aes(x = logComp1, y = logComp3, color = withinSubset, shape = withinSubset)) + geom_point( size = ptSize) + theme_bw() + 
 # layer(data = pcaFactors[pcaFactors$comp13set, ], mapping = aes(x = logComp1, y = logComp3, hjust = 0.9,vjust = -1,label = lncRnaName), geom = "text", size = 4) +
  layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp1, y = logComp3), geom = "point",size = I(ptSize)) +
 # layer(data = pcaFactors[pcaFactors$withinSubset == "true", ], mapping = aes(x = logComp1, y = logComp3, hjust = 0.9,vjust = -1,label = lncRnaName), geom = "text", size = 4) +
  ggtitle(PCA.title)+
  xlab("log(First Component)")+
  ylab("log(Third Component)")
if (labelsSet == TRUE){
  g2 = g2 + layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp1, y = logComp3, hjust = 0.9,vjust=-1,label=sub(x=lncRnaName,pattern="ENSG00000",replacement="")),geom = "text",size = 4) +    
    layer(data = pcaFactors[pcaFactors$comp12set ,], mapping = aes(x = logComp1, y = logComp3, hjust = 0.9,vjust=-1,label=lncRnaName),geom = "text",size = 4)
}

g2;ggsave(file=oneThree,height=8,width=8)

g3 <- ggplot(pcaFactors, aes(x = logComp2, y = logComp3, color = withinSubset,shape = withinSubset)) + geom_point( size = ptSize) + theme_bw() + 
  #layer(data = pcaFactors[pcaFactors$comp23set, ], mapping = aes(x = logComp2, y = logComp3, hjust = 0.9, vjust=-1, label = lncRnaName), geom = "text", size = 4)+
   layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp2, y = logComp3), geom = "point", size = I(ptSize))+
#  layer(data = pcaFactors[pcaFactors$withinSubset == "true", ], mapping = aes(x = logComp2, y = logComp3, hjust = 0.9, vjust=-1, label = lncRnaName), geom = "text", size = 4)+
  ggtitle(PCA.title)+
  xlab("log(Second Component)")+
  ylab("log(Third Component)")

if (labelsSet == TRUE){
  g3 = g3 + layer(data = pcaFactors[pcaFactors$withinSubset == "true" ,], mapping = aes(x = logComp2, y = logComp3, hjust = 0.9,vjust=-1,label=sub(x=lncRnaName,pattern="ENSG00000",replacement="")),geom = "text",size = 4) +    
    layer(data = pcaFactors[pcaFactors$comp12set ,], mapping = aes(x = logComp2, y = logComp3, hjust = 0.9,vjust=-1,label=lncRnaName),geom = "text",size = 4)
}
g3;ggsave(file=twoThree, height=8, width=8)

}


plotPcaCompToFileLabelAndWO <- function(pcaFactors,baseDir,fileNameMsg,titleMsg,comp1range=c(-Inf,Inf),comp2range=c(-Inf,Inf),comp3range=c(-Inf,Inf)){
plotPcaCompToFile(pcaFactors=pcaFactors, baseDir=baseDir, fileNameMsg=fileNameMsg, titleMsg=titleMsg,
                  comp1range=comp1range,comp2range=comp1range,comp3range=comp1range )
plotPcaCompToFile(pcaFactors=pcaFactors, baseDir=baseDir, fileNameMsg=paste(fileNameMsg,"label-",sep=""), titleMsg=titleMsg,
                  comp1range=comp1range,comp2range=comp2range,comp3range=comp3range ,labelsSet=TRUE,ptSize=I(4))
}

plotPcaPredictStats <- function(lncDf,pca,exprCols,baseDir,fileNameMsg,titleMsg){
  
  ratioTest.df <- lncDf
  ratioTest.df[exprCols] = transformDataByPCA(df=ratioTest.df[exprCols],pca=pca)
  colnames(ratioTest.df)[exprCols] = paste0("Comp.",seq_along(exprCols))
  dfo.stats <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols)
  dfo.stats.rand <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols,rand=TRUE)
  dfo.stats[["expr"]] = "realGroups"
  dfo.stats.rand[["expr"]] = "randomGroups"
  dfo.cols = c("expr","prec","sens","predict")
  dfo.comb = rbind(dfo.stats,dfo.stats.rand)
  
  
  ggplot(dfo.comb,aes(y=prec,x=sens,color=expr))+geom_point(size=1)+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth(aes(color=expr))+
         ggtitle(paste("PR Curve of radius based prediction on pca\n",titleMsg,sep=""))
  ggsave(file=paste(baseDir,fileNameMsg,"prCurve.pdf",sep=""))
  
  ggplot(dfo.comb,aes(x=FPR,y=TPR,color=expr))+geom_line()+geom_point()+theme_bw()+geom_abline(slope=1)+
    ggtitle(paste("ROC Curve of radius based prediction on pca\n",titleMsg,sep=""))
  ggsave(file=paste(baseDir,fileNameMsg,"RocCurve.pdf",sep=""))
  
  
  dfo.stats.melt <- melt(dfo.stats[c("predict","sens","prec","TPR","FPR","FNR","TNR","Fscore","accuracy","errorRate")],id.var="predict")
  dfo.stats.rand.melt <- melt(dfo.stats.rand[c("predict","sens","prec","TPR","FPR","FNR","TNR","Fscore","accuracy","errorRate")],id.var="predict")
  
}


#plot a hexbin representation of lncRNA expression components
runStatHexBin <- function(pcaFactors,baseDir,fileNameMsg,titleMsg){
  oneTwo = paste(baseDir,"comps-1-2-",fileNameMsg,"-density.pdf",sep="")  
  oneThree = paste(baseDir,"comps-1-3-",fileNameMsg,"-density.pdf",sep="")  
  twoThree = paste(baseDir,"comps-2-3-",fileNameMsg,"-density.pdf",sep="")  
  
  
  ggplot(pcaFactors,aes(x=log(Comp.1),y=log(Comp.2))) + stat_binhex()+theme_bw()+ggtitle(titleMsg)
  ggsave(file=oneTwo,height=8,width=8)
  
  ggplot(pcaFactors,aes(x=log(Comp.1),y=log(Comp.3))) + stat_binhex()+theme_bw()+ggtitle(titleMsg)
  ggsave(file=oneThree,height=8,width=8)
  
  ggplot(pcaFactors,aes(x=log(Comp.2),y=log(Comp.3))) + stat_binhex()+theme_bw()+ggtitle(titleMsg)
  ggsave(file=twoThree,height=8,width=8)
}



pcaAnalysisRemoveOutliers <- function(lncDf,exprCols,normalFun=FALSE,outDir,fileBase,foundColword,COR=FALSE){
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner <<- function(x)paste(paste("lncRNA subsection = ",foundColword),x,sep="\n")
  expr.cols <- exprCols
  
  
  lnc.pca             <- princomp(lncDf[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-loadings-all.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-scree-all.pdf"),
                          outfileECDF=makeOutFile("allLnc/lncRNA-func-scree-ECDF-all.pdf"),titleMsg="all lncRNAs")
    
  lnc.pca.factors = normalizePcaFactors(pca=lnc.pca,lncDataFrame=lncDf)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg="All Lnc RNAs",
                    comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5))
  plotPcaCompToFile(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-label-"),fileNameMsg="allLncs",titleMsg="All Lnc RNAs",
                    comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5),labelsSet=TRUE,ptSize=I(4))
  runStatHexBin(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg="PCA of lncRNA expr :\nAll Lnc RNAs")
  
  
  lncDf.reduced.1 <- lncDf[which(!lncDf$lncRnaName %in% c("MALAT1","H19","RP11-255B23.3")),]
  write.csv(file=makeOutFile("removeOut1/missingLnc.csv"),x=data.frame(missing=c("MALAT1","H19","RP11-255B23.3")))
  lnc.pca.r1             <- princomp(lncDf.reduced.1[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-loadings-removeOut1.png"),titleMsg="all lnc RNAs")
     plotScreeToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-scree-removeOut1.pdf"),
                                    outfileECDF=makeOutFile("removeOut1/lncRNA-func-scree-ECDF-removeOut1.pdf"),titleMsg="no Hg19/MALAT/RP11")
  
  lnc.pca.factors.r1 = normalizePcaFactors(pca=lnc.pca.r1,lncDataFrame=lncDf.reduced.1)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-label-"), fileNameMsg="removeOut1", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) ,labelsSet=TRUE,ptSize=I(4))
      runStatHexBin(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg="PCA of lncRNA expr :\nno Hg19/MALAT/RP11")
  
  r2.vec = c("ENSG00000235162","ENSG00000228474","ENSG00000175061","DANCR","SNGG1","ZFAS1","ENSG00000249790",
           "ENSG00000256329","ENSG00000235162","ENSG00000228474","ENSG00000249532","ENSG00000249502","MALAT1","H19","RP11-255B23.3")
  write.csv(file=makeOutFile("removeOut2/missingLnc.csv"),x=data.frame(missing=r2.vec))
  lncDf.reduced.2 <- lncDf[which(!lncDf$lncRnaName %in% r2.vec),]
  lnc.pca.r2             <- princomp(lncDf.reduced.2[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-loadings-removeOut2.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-scree-removeOut2.pdf"),
                              outfileECDF=makeOutFile("removeOut2/lncRNA-func-scree-ECDF-removeOut2.pdf"),titleMsg="group 1 & 2 of outliers removed")
  
  lnc.pca.factors.r2 = normalizePcaFactors(pca=lnc.pca.r2,lncDataFrame=lncDf.reduced.2)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11 + group 2  ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-label-"), fileNameMsg="removeOut2", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11 + group 2  ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3),labelsSet=TRUE,ptSize=I(4) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg="PCA of lncRNA expr :\nno Hg19/MALAT/RP11 + group 2")
  
  
  
  lnc.pca.factors.r2.full <- normalizePcaFactorsFromLoading(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg="Lnc RNAs: All \nLoading: no Hg19/MALAT/RP11 + group 2\n ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-label-"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg="Lnc RNAs: All \nLoading: no Hg19/MALAT/RP11 + group 2  ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3),labelsSet=TRUE,ptSize=I(4) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg="PCA of lncRNA expr :\nnLoading: o Hg19/MALAT/RP11 + group 2")
  
 
  plotPcaPredictStats(lncDf=lncDf,pca=lnc.pca.r2,exprCols=exprCols,
                    baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-PCA"),fileNameMsg="",titleMsg="removeOut2 Vectors, all lncs included")
 
  ratioTest.df <- lncDfs
  ratioTest.df[exprCols] = transformDataByPCA(df=ratioTest.df[exprCols],pca=lnc.pca.r2)
  colnames(ratioTest.df)[exprCols] = paste0("Comp.",seq_along(exprCols))
  g <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols)
  g1 <- g + ggtitle("PR Curve of radius based prediction all lncRNAs\nPCA from removeOut2 group")
  ggsave(file=makeOutFile("allLnc-removeOut2-Loadings/ratioTest.pdf"))
  
  plotPcaPredictStats(lncDf=lncDf,pca=lnc.pca.r2,exprCols=exprCols,
                      baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-PCA"),fileNameMsg="",titleMsg="removeOut2 Vectors, all lncs included")
  


}

pcaAnalysisRemoveOutliersAllCols  <- function(lncDf,exprCols,normalFun=FALSE,outDir,fileBase,foundColword,COR=FALSE){
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner <<- function(x)paste(paste("lncRNA subsection = ",foundColword),x,sep="\n")
  expr.cols <- exprCols
  
  
  lnc.pca             <- princomp(lncDf[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-loadings-all.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-scree-all.pdf"),
                  outfileECDF=makeOutFile("allLnc/lncRNA-func-scree-ECDF-all.pdf"),titleMsg="all lncRNAs")
  
  lnc.pca.factors = normalizePcaFactors(pca=lnc.pca,lncDataFrame=lncDf)
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg="All Lnc RNAs",
                    comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5))

  runStatHexBin(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg="PCA of lncRNA expr :\nAll Lnc RNAs")
  
  
  lncDf.reduced.1 <- lncDf[which(!lncDf$lncRnaName %in% c("MALAT1","H19","RP11-255B23.3")),]
  write.csv(file=makeOutFile("removeOut1/missingLnc.csv"),x=data.frame(missing=c("MALAT1","H19","RP11-255B23.3")))
  lnc.pca.r1             <- princomp(lncDf.reduced.1[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-loadings-removeOut1.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-scree-removeOut1.pdf"),
                  outfileECDF=makeOutFile("removeOut1/lncRNA-func-scree-ECDF-removeOut1.pdf"),titleMsg="no Hg19/MALAT/RP11")
  
 
  lnc.pca.factors.r1 = normalizePcaFactors(pca=lnc.pca.r1,lncDataFrame=lncDf.reduced.1)
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
   runStatHexBin(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg="PCA of lncRNA expr :\nno Hg19/MALAT/RP11")
  
  r2.vec = c("ENSG00000235162","ENSG00000228474","ENSG00000175061","DANCR","SNGG1","ZFAS1","ENSG00000249790",
             "ENSG00000256329","ENSG00000235162","ENSG00000228474","ENSG00000249532","ENSG00000249502","MALAT1","H19","RP11-255B23.3")
  write.csv(file=makeOutFile("removeOut2/missingLnc.csv"),x=data.frame(missing=r2.vec))
  lncDf.reduced.2 <- lncDf[which(!lncDf$lncRnaName %in% r2.vec),]
  lnc.pca.r2             <- princomp(lncDf.reduced.2[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-loadings-removeOut2.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-scree-removeOut2.pdf"),
                  outfileECDF=makeOutFile("removeOut2/lncRNA-func-scree-ECDF-removeOut2.pdf"),titleMsg="group 1 & 2 of outliers removed")
  
  lnc.pca.factors.r2 = normalizePcaFactors(pca=lnc.pca.r2,lncDataFrame=lncDf.reduced.2)
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func"), fileNameMsg="removeOut2", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11 + group 2  ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg="PCA of lncRNA expr :\nno Hg19/MALAT/RP11 + group 2")
  
  lnc.pca.factors.r2.full <- normalizePcaFactorsFromLoading(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols)
  exprortAsTable(df=lnc.pca.factors.r2.full,file=makeOutFile("removeOut2-loadings/pca.tab"))  
  
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg="Lnc RNAs: All \nLoading: no Hg19/MALAT/RP11 + group 2\n ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
   runStatHexBin(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg="PCA of lncRNA expr :\nnLoading: o Hg19/MALAT/RP11 + group 2")
  
  ratioTest.df <- lncDf
  ratioTest.df[exprCols] = transformDataByPCA(df=ratioTest.df[exprCols],pca=lnc.pca.r2)
  colnames(ratioTest.df)[exprCols] = paste0("Comp.",seq_along(exprCols))
  g <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols)
  g1 <- g + ggtitle("PR Curve of radius based prediction all lncRNAs\nPCA from removeOut2 group")
  ggsave(file=makeOutFile("allLnc-removeOut2-Loadings/ratioTest.pdf"))
}
pcaAnalysisRemoveOutliersGeneral  <- function(lncDf,exprCols,titleMsg,normalFun=FALSE,outDir,fileBase,foundColword,COR=FALSE,rOneVec=c("MALAT1","H19","RP11-255B23.3"),
                                              rTwoVec=c("ENSG00000235162","ENSG00000228474","ENSG00000175061","DANCR","SNGG1","ZFAS1","ENSG00000249790",
                                                        "ENSG00000256329","ENSG00000235162","ENSG00000228474","ENSG00000249532","ENSG00000249502","MALAT1","H19","RP11-255B23.3") ){
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner <<- function(x)paste(paste("lncRNA subsection = ",foundColword),x,sep="\n")
  expr.cols <- exprCols
  
  dirs = c(makeOutFile("allLnc/"),makeOutFile("removeOut1/"),makeOutFile("removeOut2/"),makeOutFile("allLnc-removeOut2-Loadings/"))
  for (newDir in dirs){
    if (!file.exists(newDir)){
      dir.create(newDir)
    }
  }
  
  
  
  lnc.pca             <- princomp(lncDf[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-loadings-all.png"),titleMsg=paste("all lnc RNAs",titleMsg))
  plotScreeToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-scree-all.pdf"),
                  outfileECDF=makeOutFile("allLnc/lncRNA-func-scree-ECDF-all.pdf"),titleMsg=paste("all lnc RNAs",titleMsg))
  

  lnc.pca.factors = normalizePcaFactors(pca=lnc.pca,lncDataFrame=lncDf)
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg=paste("all lnc RNAs",titleMsg),
                              comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5))
  
  runStatHexBin(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg=paste("PCA of lncRNA expr :\n",titleMsg,sep=""))
  
  
  lncDf.reduced.1 <- lncDf[which(!lncDf$lncRnaName %in% c("MALAT1","H19","RP11-255B23.3")),]
  write.csv(file=makeOutFile("removeOut1/missingLnc.csv"),x=data.frame(missing=c("MALAT1","H19","RP11-255B23.3")))
  lnc.pca.r1             <- princomp(lncDf.reduced.1[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-loadings-removeOut1.png"),titleMsg=paste("all lnc RNAs",titleMsg))
  plotScreeToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-scree-removeOut1.pdf"),
                  outfileECDF=makeOutFile("removeOut1/lncRNA-func-scree-ECDF-removeOut1.pdf"),titleMsg=paste("no Hg19/MALAT/RP11",titleMsg))
  
  
  lnc.pca.factors.r1 = normalizePcaFactors(pca=lnc.pca.r1,lncDataFrame=lncDf.reduced.1)
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg=paste("Lnc RNAs: no Hg19/MALAT/RP11",titleMsg),
                              comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg=paste("PCA of lncRNA expr :\nno Hg19/MALAT/RP11",titleMsg))
  
  r2.vec = rTwoVec
  write.csv(file=makeOutFile("removeOut2/missingLnc.csv"),x=data.frame(missing=r2.vec))
  lncDf.reduced.2 <- lncDf[which(!lncDf$lncRnaName %in% r2.vec),]
  lnc.pca.r2             <- princomp(lncDf.reduced.2[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-loadings-removeOut2.png"),titleMsg=paste("all lnc RNAs",titleMsg))
  plotScreeToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-scree-removeOut2.pdf"),
                  outfileECDF=makeOutFile("removeOut2/lncRNA-func-scree-ECDF-removeOut2.pdf"),titleMsg=paste("group 1 & 2 of outliers removed",titleMsg))
  
  lnc.pca.factors.r2 = normalizePcaFactors(pca=lnc.pca.r2,lncDataFrame=lncDf.reduced.2)
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func"), fileNameMsg="removeOut2", titleMsg=paste("Lnc RNAs: no Hg19/MALAT/RP11 + group 2  \n",titleMsg),
                              comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg=paste("PCA of lncRNA expr :\nno Hg19/MALAT/RP11 + group 2\n",titleMsg))
  
  lnc.pca.factors.r2.full <- normalizePcaFactorsFromLoading(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols)
  
  
  #normalizePcaFactorsFromLoadingForRatioTestg(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols)
  exportAsTable(df= normalizePcaFactorsFromLoadingForRatioTest(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols),
                file=makeOutFile("pca.tab"))  
  
  print("made it to removeOut2-Loadings...")
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg=paste("Lnc RNAs: All \nLoading: no Hg19/MALAT/RP11 + group 2\n ",titleMsg),
                              comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg=paste("PCA of lncRNA expr :\nnLoading: o Hg19/MALAT/RP11 + group 2",titleMsg) )
  
  print("made it to the radius test")
#  ratioTest.df <- lncDf
#  ratioTest.df[exprCols] = transformDataByPCA(df=ratioTest.df[exprCols],pca=lnc.pca.r2)
#  colnames(ratioTest.df)[exprCols] = paste0("Comp.",seq_along(exprCols))

 # g <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols) + ggtitle("PR Curve of radius based prediction all lncRNAs\nPCA from removeOut2 group")

 # g 
 # ggsave(file=makeOutFile("allLnc-removeOut2-Loadings/ratioTest.pdf"))
  
 #plotPcaPredictStats(lncDf=lncDf,pca=lnc.pca.r2,exprCols=exprCols,
#                      baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-PCA"),fileNameMsg="",titleMsg="removeOut2 Vectors, all lncs included")
  
  
}

pcaAnalysisRemoveOutliersRobust <- function(lncDf,exprCols,normalFun=FALSE,outDir,fileBase,foundColword,COR=FALSE){
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner <<- function(x)paste(paste("lncRNA subsection = ",foundColword),x,sep="\n")
  expr.cols <- exprCols
  
  
  lnc.pca             <- princomp(lncDf[exprCols],covmat=MASS::cov.rob(lncDf[exprCols]))
  plotLoadingsToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-loadings-all.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-scree-all.pdf"),
                  outfileECDF=makeOutFile("allLnc/lncRNA-func-scree-ECDF-all.pdf"),titleMsg="all lncRNAs")
  
  lnc.pca.factors = normalizePcaFactors(pca=lnc.pca,lncDataFrame=lncDf)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg="All Lnc RNAs",
                    comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5))
  plotPcaCompToFile(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-label-"),fileNameMsg="allLncs",titleMsg="All Lnc RNAs",
                    comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5),labelsSet=TRUE,ptSize=I(4))
  runStatHexBin(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg="PCA of lncRNA expr :\nAll Lnc RNAs")
  
  
  lncDf.reduced.1 <- lncDf[which(!lncDf$lncRnaName %in% c("MALAT1","H19","RP11-255B23.3")),]
  write.csv(file=makeOutFile("removeOut1/missingLnc.csv"),x=data.frame(missing=c("MALAT1","H19","RP11-255B23.3")))
  lnc.pca.r1             <- princomp(lncDf.reduced.1[exprCols],covmat=MASS::cov.rob(lncDf.reduced.1[exprCols]))
  plotLoadingsToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-loadings-removeOut1.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-scree-removeOut1.pdf"),
                  outfileECDF=makeOutFile("removeOut1/lncRNA-func-scree-ECDF-removeOut1.pdf"),titleMsg="no Hg19/MALAT/RP11")
  
  lnc.pca.factors.r1 = normalizePcaFactors(pca=lnc.pca.r1,lncDataFrame=lncDf.reduced.1)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-label-"), fileNameMsg="removeOut1", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) ,labelsSet=TRUE,ptSize=I(4))
  runStatHexBin(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg="PCA of lncRNA expr :\nno Hg19/MALAT/RP11")
  
  r2.vec = c("ENSG00000235162","ENSG00000228474","ENSG00000175061","DANCR","SNGG1","ZFAS1","ENSG00000249790",
             "ENSG00000256329","ENSG00000235162","ENSG00000228474","ENSG00000249532","ENSG00000249502","MALAT1","H19","RP11-255B23.3")
  write.csv(file=makeOutFile("removeOut2/missingLnc.csv"),x=data.frame(missing=r2.vec))
  lncDf.reduced.2 <- lncDf[which(!lncDf$lncRnaName %in% r2.vec),]
  lnc.pca.r2             <- princomp(lncDf.reduced.2[exprCols],covmat=MASS::cov.rob(lncDf.reduced.2[exprCols]))
  plotLoadingsToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-loadings-removeOut2.png"),titleMsg="all lnc RNAs")
  plotScreeToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-scree-removeOut2.pdf"),
                  outfileECDF=makeOutFile("removeOut2/lncRNA-func-scree-ECDF-removeOut2.pdf"),titleMsg="group 1 & 2 of outliers removed")
  
  lnc.pca.factors.r2 = normalizePcaFactors(pca=lnc.pca.r2,lncDataFrame=lncDf.reduced.2)
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11 + group 2  ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  plotPcaCompToFile(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-label-"), fileNameMsg="removeOut2", titleMsg="Lnc RNAs: no Hg19/MALAT/RP11 + group 2  ",
                    comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3),labelsSet=TRUE,ptSize=I(4) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg="PCA of lncRNA expr :\nno Hg19/MALAT/RP11 + group 2")
  
  ratioTest.df <- lncDf
  ratioTest.df[exprCols] = transformDataByPCA(df=ratioTest.df[exprCols],pca=lnc.pca.r2)
  colnames(ratioTest.df)[exprCols] = paste0("Comp.",seq_along(exprCols))
  g <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols)
  g1 <- g + ggtitle("PR Curve of radius based prediction all lncRNAs\nPCA from removeOut2 group")
  ggsave(file=makeOutFile("removeOut2/ratioTest.pdf"))
}

findTopOutliersFromPCA <- function(pca,N,exprCols){
  #pca should be normalized...
  pcaCols <- seq_along(exprCols)
  pca.center <- suppressMessages(melt(ddply(pca[pcaCols],.(),numcolwise(mean)))[["value"]])
  pca[["dist"]] <- apply(pca[pcaCols],1,function(x) sqrt(sum((x - pca.center)^2)))
  pca.distOrder <- pca[order(max(pca$dist) - pca$dist),]
  pca.distOrder[1:N,"lncRnaName"]
}



pcaAnalysisRemoveOutliersSelectOutliers  <- function(lncDf,exprCols,titleMsg,normalFun=FALSE,outDir,fileBase,foundColword,COR=FALSE ){
  printReport <- function(report){ if(verbose == TRUE){print(report)}}
  makeOutFile <- function(x){outfile<-paste(paste(outDir,fileBase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & fileBase
  titleWithBanner <<- function(x)paste(paste("lncRNA subsection = ",foundColword),x,sep="\n")
  expr.cols <- exprCols
  
  dirs = c(makeOutFile("allLnc/"),makeOutFile("removeOut1/"),makeOutFile("removeOut2/"),makeOutFile("allLnc-removeOut2-Loadings/"))
  for (newDir in dirs){
    if (!file.exists(newDir)){
      dir.create(newDir)
    }
  }
  
  
  
  ### PCA FULL
  lnc.pca             <- princomp(lncDf[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-loadings-all.png"),titleMsg=paste("all lnc RNAs",titleMsg))
  plotScreeToFile(pca=lnc.pca,outfile=makeOutFile("allLnc/lncRNA-func-scree-all.pdf"),
                  outfileECDF=makeOutFile("allLnc/lncRNA-func-scree-ECDF-all.pdf"),titleMsg=paste("all lnc RNAs",titleMsg))
  
    #normalize for plotting...
  lnc.pca.factors = normalizePcaFactors(pca=lnc.pca,lncDataFrame=lncDf)
 
  saveToFilePlotDistAway(lnc.pca.factors,seq_along(exprCols),makeOutFile("allLnc/"),"PR curve on PCA - all lncRNA")
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg=paste("all lnc RNAs",titleMsg),
                              comp1range=c(0,4),comp2range=c(2,10),comp3range=c(0,3.5))
  
  runStatHexBin(pcaFactors=lnc.pca.factors,baseDir=makeOutFile("allLnc/lncRNA-func-"),fileNameMsg="allLncs",titleMsg=paste("PCA of lncRNA expr :\n",titleMsg,sep=""))
  
  
  ### PCA Remove group 1
  ## need to find top N outlier
  r1.vec =  findTopOutliersFromPCA(lnc.pca.factors,10,expr.cols)
  lncDf.reduced.1 <- lncDf[which(!lncDf$lncRnaName %in% r1.vec),]
  write.csv(file=makeOutFile("removeOut1/missingLnc.csv"),x=data.frame(missing=r1.vec))
  lnc.pca.r1             <- princomp(lncDf.reduced.1[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-loadings-removeOut1.png"),titleMsg=paste("all lnc RNAs",titleMsg))
  plotScreeToFile(pca=lnc.pca.r1, outfile=makeOutFile("removeOut1/lncRNA-func-scree-removeOut1.pdf"),
                  outfileECDF=makeOutFile("removeOut1/lncRNA-func-scree-ECDF-removeOut1.pdf"),titleMsg=paste("no r1 group",titleMsg))
  
     #normalize for plotting...
  lnc.pca.factors.r1 = normalizePcaFactors(pca=lnc.pca.r1,lncDataFrame=lncDf.reduced.1)
  saveToFilePlotDistAway(lnc.pca.factors,seq_along(exprCols),makeOutFile("removeOut1/"),"PR curve on PCA - remove group 1")

  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg=paste("Lnc RNAs: no r1 group",titleMsg),
                              comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r1, baseDir=makeOutFile("removeOut1/lncRNA-func-"), fileNameMsg="removeOut1", titleMsg=paste("PCA of lncRNA expr :\nno  r1 group",titleMsg))
  
  
  ### PCA Remove group 2
  ## need to find top N outlier
  r2.vec =  c(findTopOutliersFromPCA(lnc.pca.factors.r1,10,expr.cols),r1.vec)
  write.csv(file=makeOutFile("removeOut2/missingLnc.csv"),x=data.frame(missing=r2.vec))
  lncDf.reduced.2 <- lncDf[which(!lncDf$lncRnaName %in% r2.vec),]
  lnc.pca.r2             <- princomp(lncDf.reduced.2[exprCols],cor=COR)
  plotLoadingsToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-loadings-removeOut2.png"),titleMsg=paste("all lnc RNAs",titleMsg))
  plotScreeToFile(pca=lnc.pca.r2, outfile=makeOutFile("removeOut2/lncRNA-func-scree-removeOut2.pdf"),
                  outfileECDF=makeOutFile("removeOut2/lncRNA-func-scree-ECDF-removeOut2.pdf"),titleMsg=paste("group 1 & 2 of outliers removed",titleMsg))
  
  lnc.pca.factors.r2 = normalizePcaFactors(pca=lnc.pca.r2,lncDataFrame=lncDf.reduced.2)
  saveToFilePlotDistAway(lnc.pca.factors.r2,seq_along(exprCols),makeOutFile("removeOut2/"),"PR curve on PCA - remove group 2")
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func"), fileNameMsg="removeOut2", titleMsg=paste("Lnc RNAs: no r1 group + r2 group   \n",titleMsg),
                              comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2, baseDir=makeOutFile("removeOut2/lncRNA-func-"), fileNameMsg="removeOut2", titleMsg=paste("PCA of lncRNA expr :\nno no r1 group + r2 group \n",titleMsg))
  
  
  
  ### PCA all lncRNA, w/ loadings from remove2 grouping...
  print("making r2 loading on all lncRNA")
  lnc.pca.factors.r2.full <- normalizePcaFactorsFromLoading(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols)
  doubleCols = as.vector(which(sapply(colnames(lnc.pca.factors.r2.full),function(x)(typeof(lnc.pca.factors.r2.full[1,x]) == "double"))))
  saveToFilePlotDistAway(lnc.pca.factors.r2.full,doubleCols,makeOutFile("removeOut2/"),"PR curve on PCA - remove group 2")
  
  
  #normalizePcaFactorsFromLoadingForRatioTestg(pca=lnc.pca.r2,lncDataFrame=lncDf,exprCols=exprCols)
  exportAsTable(df= lnc.pca.factors.r2.full,
                file=makeOutFile("allLnc-removeOut2-Loadings/pca.tab"))  
  
  print("made it to removeOut2-Loadings...")
  plotPcaCompToFileLabelAndWO(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg=paste("Lnc RNAs: All \nLoading: no Hg19/MALAT/RP11 + group 2\n ",titleMsg),
                              comp1range=c(1.4,Inf),comp2range=c(0,3.5),comp3range=c(0,3) )
  runStatHexBin(pcaFactors=lnc.pca.factors.r2.full, baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-"), fileNameMsg="allLnc-removeOut2-Loadings", titleMsg=paste("PCA of lncRNA expr :\nnLoading: o Hg19/MALAT/RP11 + group 2",titleMsg) )
  
  print("made it to the radius test")
  #  ratioTest.df <- lncDf
  #  ratioTest.df[exprCols] = transformDataByPCA(df=ratioTest.df[exprCols],pca=lnc.pca.r2)
  #  colnames(ratioTest.df)[exprCols] = paste0("Comp.",seq_along(exprCols))
  
  # g <- plotDistanceAwayRatio(df=ratioTest.df,exprCols=exprCols) + ggtitle("PR Curve of radius based prediction all lncRNAs\nPCA from removeOut2 group")
  
  # g 
  # ggsave(file=makeOutFile("allLnc-removeOut2-Loadings/ratioTest.pdf"))
  
  # plotPcaPredictStats(lncDf=lncDf,pca=lnc.pca.r2,exprCols=exprCols,
  #                     baseDir=makeOutFile("allLnc-removeOut2-Loadings/lncRNA-func-PCA"),fileNameMsg="",titleMsg="removeOut2 Vectors, all lncs included")
  
  
}


main <- function(){
  
combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
lnc.out.dir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs"
lnc.out.dir.robust = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs-robust"

lnc.colIndex = 2:33
df.1=readInTable(lnc.in.file)
df.1$gene_id_short <- sapply(df.1$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
func.df <- getEnslist()

func.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs"

print("start functional List")

func.df <- within( func.df, {
  gene_id_short = ensembl_gene_id
  gene_id = gene_id_short
  lncRnaName = external_gene_id})


f.df <- getAnnotLncDf(lncDf=df.1,annotDf=func.df,exprCol=2:33,annotColName="lncRnaName")
f.df[["lncRnaName"]] = ifelse(f.df$lncRnaName == "notFound",f.df$gene_id_short,f.df$lncRnaName)
pcaAnalysisRemoveOutliersGeneral(lncDf=f.df,exprCols=2:33,foundColword="functional-LncRNA",fileBase="",outDir=func.outdir)
#source("/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/pcaAnalysisWOoutlier.R");pcaAnalysisRemoveOutliers(lncDf=f.df,exprCols=2:33,foundColword="functional-LncRNA",fileBase="",outDir=func.outdir)


#lnc.out.dir.robust = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs-robust"
#m = mean(f.df[[x]])
#s = sd(f.df[[x]])
#mi = (min(f.df[[x]]) - m) / s
#robust.df <- as.data.frame(sapply(colnames(f.df[,2:33]),function(x){((f.df[[x]] - mean(f.df[[x]])) / sd(f.df[[x]])) + (min(f.df[[x]]) - mean(f.df[[x]])) / sd(f.df[[x]]) + 1  }))
#pcaAnalysisRemoveOutliersRobust(lncDf=robust.df,exprCols=1:32,foundColword="functional-LncRNA",fileBase="PCA-",outDir=lnc.out.dir.robust)

lnc.out.dir.allCols = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs-allColsNorm"
func.outdir        = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs"
lpa.outdir        = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs-longPolyA"
lnpa.outdir        = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs-longNonPolyA"

num.cols = as.vector(which(unlist(sapply(colnames(f.df), function(x)(typeof(f.df[[x]]) == "double")))))
f.df$nTimesRest  = ifelse(f.df$nTimesRest == Inf, 0,f.df$nTimesRest)
norm <- as.data.frame(sapply(colnames(f.df[,num.cols]),function(x){(f.df[[x]] - mean(f.df[[x]])) /sd(f.df[[x]])}))
norm.df <- cbind(norm,f.df[,which(!colnames(f.df) %in% colnames(norm))])

pcaAnalysisRemoveOutliersGeneral(lncDf=f.df,exprCols=2:33,foundColword="functional-LncRNA",fileBase="",outDir=func.outdir,titleMsg="lpa + lnpa rna seq data")
pcaAnalysisRemoveOutliersGeneral(lncDf=f.df,exprCols=c(num.cols)[-34], # remove sumExpr, as it is equiv to average expression...(linearly dependant)
                                 foundColword="functional-LncRNA",fileBase="",outDir=lnc.out.dir.allCols,titleMsg="lpa + lnpa rna seq data+\nall other cols")

pcaAnalysisRemoveOutliersGeneral(lncDf=f.df,exprCols=grep("longNonPolyA$",colnames(f.df)),foundColword="functional-LncRNA",fileBase="",outDir=lpa.outdir,titleMsg="longNonPolyA rna-seq")
pcaAnalysisRemoveOutliersGeneral(lncDf=f.df,exprCols=grep("longPolyA$",colnames(f.df)),foundColword="functional-LncRNA",fileBase="",outDir=lpa.outdir,titleMsg="longPolyA rna-seq")
}

#  grep("longNonPolyA$",colnames(f.df))
#

