#The purpose of this script is to run: lncRNA comparison, PCA, logistic Regression, and eigenRank 
# the output will be directed to one directory
#
#
#
home <- Sys.getenv("HOME")

getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)






source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/eigenRank_08272013.R",sep=""))
source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/rnaSeq_qc.R",sep=""))
source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/logisticRegPcaExprData.R",sep=""))
#source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R",sep=""))
library(class)
library(MASS)
library(LiblineaR)
library(ROCR) # http://cran.r-project.org/web/packages/ROCR/index.html
library(glmnet) # http://cran.r-project.org/web/packages/glmnet/glmnet.pdf

#LOG Reg 



runLogRegIter <- function(lncDf=local.df,cols=exprLogReg,outdir,msg,trials=10){
  if(sum(local.df$label) < 3){
    return(0)
  }
  
  
  if (!file.exists(outdir)){
    dir.create(outdir,recursive=TRUE)
  }
  
  ratio.df <- expand.grid(mainEffects = c(TRUE),reg = c(FALSE),ratio=seq(0.7),run=1:10)
  ratio.df$predict <- NA

  
  print(outdir)
  for(i in 1:dim(ratio.df)[1]){
    print(i)
    
    line <- ratio.df[i,]
    output <- trainAndTestLogReg(lncDf=lncDf,cols=cols, ratio = line$ratio,
                                 mainEffects=line$mainEffects,reg=line$reg)
    
    
    ratio.df$predict[i] <- output$R
    ratio.df$TPR[i] <- output$TPR
    ratio.df$TP[i] <- output$TP
    ratio.df$FP[i] <- output$FP
    ratio.df$FN[i] <- output$FN
    ratio.df$TN[i] <- output$TN
    ratio.df$AUC[i] <- output$AUC
  }
  exportAsTable(ratio.df,file= paste(outdir,"ratioExpr.tab",sep="/"))
  #lambda exper
  lambda.df <- as.data.frame(expand.grid(mainEffects = c(TRUE),ratio=0.7,reg = c(TRUE),lambda=c(0,10 ^ (seq(-2,8))),run=1:10))
  lambda.df$predict <- NA
  for(i in 1:dim(lambda.df)[1]){
    print(i)
    line <- lambda.df[i,]
    output <- trainAndTestLogReg(lncDf = lncDf, cols=cols, ratio = line$ratio, 
                                 mainEffects = line$mainEffects, reg = line$reg, lambda = line$lambda)
    lambda.df$predict[i] <- output$R
    lambda.df$TPR[i] <- output$TPR
    
    lambda.df$TP[i] <- output$TP
    lambda.df$FP[i] <- output$FP
    lambda.df$FN[i] <- output$FN
    lambda.df$TN[i] <- output$TN
    lambda.df$AUC[i] <- output$AUC
    
  }
  exportAsTable(lambda.df, file=paste(outdir,"lambdaExpr.tab",sep="/"))
  analyzeLogRegTest(outdir, indir=outdir,msg=msg)
}

plotLncRNAComparisons <- function(lncDf=local.df,outdir=outdir,  cols=expr.cols ,titleMsg="",
                                  foundColword="",filebase=""){
  
  
  if(!file.exists(outdir)){
    dir.create(path =outdir,recursive=TRUE)
    
  }
  lncDf$withinSubset = ifelse(lncDf$label == 1, "true","false")
  n.subset = length(which(lncDf$withinSubset == "true"))
  n.total = dim(lncDf)[1]
  n.line  = paste("subset =",foundColword,"\n",n.subset,foundColword,"lncRNAs out of",n.total,"total lncRNAs\n",sep=" ")
  titleWithBanner = function(x)paste(titleMsg,x,sep="\n")
  makeOutFile <- function(x){outfile<-paste(paste(outdir,filebase,sep="/"),x,sep="");print(paste("making",outfile));outfile} # requires outDir & filebase
  
  
  lncDf <- editStatsForLncDf(expr.df = lncDf, cols = cols)
  
  annotColName <- "withinSubset"
  id.vec <- c("entropyExpr","threeTimesRest","maxExprType","tissSpec","withinSubset")
  melt.df <- melt(lncDf[c("gene_id",cols,id.vec)], 
                  measure.vars=sort(cols),
                  id.vars=c("gene_id",id.vec))
  
  
  
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
  
  
  
  
  #polar plot 
  ggplot(melt.df,aes(x=variable,y=value,group=gene_id,color=factor(withinSubset))) + 
    geom_line()+coord_polar()+
    theme_bw()+
    ggtitle(titleWithBanner("star plot of total expression over all samples"))
  ggsave(file=makeOutFile("polar.pdf"),height=16,width=16)
  
  
  top10LncRNA = unique(melt.df[order(melt.df$value,decreasing=TRUE),"gene_id"])[1:10]
  ggplot(melt.df[which(!melt.df$gene_id %in% top10LncRNA),],aes(x=variable,y=value,group=gene_id,color=factor(withinSubset))) + 
    geom_line()+coord_polar()+
    theme_bw()+
    ggtitle(titleWithBanner("star plot of total expression over all samples\nwithout top 10 highest lncRNA"))
  ggsave(file=makeOutFile("polar-removTop10.pdf"),height=16,width=16)
  
  ggplot(melt.df,aes(x=variable,y=log2(value),group=gene_id,color=factor(withinSubset))) + 
    geom_line()+coord_polar()+
    theme_bw()+
    ggtitle(titleWithBanner("star plot of total expression over all samples"))
  ggsave(file=makeOutFile("log-polar.pdf"),height=16,width=16)
  
  ggplot(lncDf, aes(x= averageExpr, y = tissSpec, color = factor(withinSubset)))+
    geom_point() +
    theme_bw() +
    ggtitle(titleWithBanner("average expression vs. tissue specificity (JSD)"))
  ggsave(file=makeOutFile("tissSpec-vs-aveExpr.pdf"),height=16,width=16)
  
  
  ggplot(lncDf, aes(x= log(averageExpr), y = tissSpec, color = factor(withinSubset),size=label))+
    geom_point() +
    theme_bw() +
    ggtitle(titleWithBanner("average expression vs. tissue specificity (JSD)"))
  ggsave(file=makeOutFile("log-tissSpec-vs-aveExpr.pdf"),height=16,width=16)
  
  foundInExpr.df = as.data.frame(table(melt.df[which(melt.df$value > 0 ), "maxExprType"]))
  colnames(foundInExpr.df) <- c("expr", "count")
  ggplot(foundInExpr.df, aes(x = expr, y = count))+ geom_bar(stat="identity")+
    theme_bw()+
    coord_flip()
  ggsave(file=makeOutFile("lncRNA-foundInEachExpr.pdf"),height=9,width=6)
  
  
  ggplot(lncDf,aes(x=tissSpec,fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
    ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
    xlab("JSD for each gene")
  ggsave(file=makeOutFile( "lncRNA-JSD-density-plot.pdf"),height=8,width=8)
  
  ggplot(lncDf,aes(x=tissSpec,fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+ facet_wrap(~withinSubset)+
    ggtitle(titleWithBanner("expression accross all samples(cell-type/pulldown)"))+
    xlab("JSD for each gene")
  ggsave(file=makeOutFile("lncRNA-JSD-density-plot.pdf"),height=8,width=8)
  
  
  # Comparison by 
  # density, freqpoly
  # log and nolog
  ggplot(stat.melt,aes(x=log(value),fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("measures across all samples(cell-type/pulldown)"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("log-density-plot.pdf"),height=8,width=8)
  
  ggplot(stat.melt,aes(x=log(value),fill=withinSubset))+geom_freqpoly() + theme_bw()+
    facet_grid(withinSubset~variable,scale="free") +
    ggtitle(titleWithBanner("measures accross all samples(cell-type/pulldown)"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("log-allcols-freqpoly-plot.pdf"),height=6,width=13)
  
  ggplot(stat.melt,aes(x=value,fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("measures accross all samples(cell-type/pulldown)"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("compare-density-plot.pdf"),height=8,width=8)
  
  ggplot(stat.melt,aes(x=value,fill=withinSubset))+geom_freqpoly() + theme_bw()+
    facet_grid(withinSubset~variable,scale="free") +
    ggtitle(titleWithBanner("measures accross all samples(cell-type/pulldown)"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("compare-freqpoly-plot.pdf"),height=6,width=13)
  
  
  
  ## Box plots
  # log and no log
  ggplot(stat.melt,aes(x=withinSubset,y=log(value),fill=withinSubset))+geom_boxplot() + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("measures across all samples(cell-type/pulldown)"))+
    xlab("value of facet label")
  #+scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-allcols-log-box-plot.pdf"),height=8,width=8)
  
  ggplot(stat.melt,aes(x=withinSubset,y=value,fill=withinSubset))+geom_boxplot() + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("measures across all samples(cell-type/pulldown)"))+
    xlab("value of facet label")
  #+scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-allcols-box-plot.pdf"),height=8,width=8)
  
  
  ## Derrien information plots
  
  ggplot(derr.melt,aes(x=log(value),fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-derrien-log-density-plot.pdf"),height=8,width=8)
  
  ggplot(derr.melt,aes(x=log(value),fill=withinSubset))+geom_freqpoly() + theme_bw()+
    facet_grid(withinSubset~variable,scale="free") +
    ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-derrien-log-freqpoly-plot.pdf"),height=6,width=16)
  
  ggplot(derr.melt,aes(x=withinSubset,y=log(value),fill=withinSubset))+geom_boxplot() + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
    xlab("value of facet label")
  #+scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-derrien-log-box-plot.pdf"),height=8,width=8)
  
  ggplot(derr.melt,aes(x=value,fill=withinSubset))+geom_density(alpha=I(0.4)) + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-derrien-allcols-density-plot.pdf"),height=8,width=8)
  
  ggplot(derr.melt,aes(x=value,fill=withinSubset))+geom_freqpoly() + theme_bw()+
    facet_grid(withinSubset~variable,scale="free") +
    ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
    xlab("value of facet label")+
    scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-derrien-allcols-freqpoly-plot.pdf"),height=6,width=16)
  
  ggplot(derr.melt,aes(x=withinSubset,y=value,fill=withinSubset))+geom_boxplot() + theme_bw()+
    facet_wrap(~variable,scale="free") +
    ggtitle(titleWithBanner("average features from Derrien 2012 datasheet"))+
    xlab("value of facet label")
  #+scale_colour_manual(values = c("true"="green","false"="red"))
  ggsave(file=makeOutFile("lncRNA-derrien-allcols-box-plot.pdf"),height=8,width=8)
  
}



#plotLncRNAComparisons <- function(lncDf=local.df,outdir=outdir,
#                                  cols=expr.cols ,titleMsg="",
#                                  foundColword="",filebase=""){


runPCA_helper <- function(lncDf=local.df,outdir=outdir,
                          cols=expr.cols, titleMsg="",
                          foundColword="",filebase=""){
  
  func.df <- getEnslist()
  func.outdir = "/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/PCA-funcLncs"
  
  print("start functional List")
  
  func.df <- within( func.df, {
    gene_id_short = ensembl_gene_id
    gene_id = gene_id_short
    lncRnaName = external_gene_id})
  
  
  f.df <- getAnnotLncDf(lncDf=df.1,annotDf=func.df,exprCol=2:33,annotColName="lncRnaName")
  f.df[["lncRnaName"]] = ifelse(f.df$lncRnaName == "notFound",f.df$gene_id_short,f.df$lncRnaName)
  
  pcaAnalysisRemoveOutliersGeneral(lncDf=lncDf,exprCols=expr.cols,foundColword="functional-LncRNA",filebase=filebase,outDir=outdir)
  
  
  
  
}




#################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

# rm(list=unlist(ls()))
# source('~/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/fullAnalysis-Jan2014.R')
# runCompleteAnalysis(); runLogRegTest()

# library(doParallel);library(foreach)
# registerDoParallel(10)
# df <-  expand.grid(c("one","two"),c("a","b","c","d"))
# one <- foreach(i = dim(df)[1]) %dopar% {Sys.sleep(runif(1));print(i[i,]);print("end") }
# Reduce(x=seq_along(one), function(x,y){cbind(x,one[[y]])})
runCompleteAnalysis <-function(outdirFull =  getFullPath("plots/fullAnalysisExperiment/"),
                               lncGeneIdShorts.vec = c("IDRlessthan0_01", "IDRlessthan0_1","IDRlessthan0_2","allLncRNA"),
                               biotypes.vec = c("antisense","lincRNA","processed_transcript","all_biotypes","remove_antisense","lincProcTrans"),
                               cols.vec = c("lpa","bothPullDowns","lnpa")){
  
  bm73.file= getFullPath("data/lnc_biotype_ens73.tab")
  bm.df <- readInTable(file=bm73.file)
  bm.df$bm_biotype <- bm.df$gene_biotype
  bm.df$gene_id_short <- bm.df$ensembl_gene_id
  
  biotypes.list = list(antisense = bm.df[which(bm.df$bm_biotype == "antisense"),"gene_id_short"],
                       processed_transcript = bm.df[which(bm.df$bm_biotype == "processed_transcript"),"gene_id_short"],
                       lincRNA = bm.df[which(bm.df$bm_biotype == "lincRNA"),"gene_id_short"],
                       all_biotypes = bm.df[["gene_id_short"]],
                       remove_antisense =  bm.df[which(bm.df$bm_biotype != "antisense"),"gene_id_short"],
                       lincProcTrans = bm.df[which(bm.df$bm_biotype == "lincRNA" | bm.df$bm_biotype == "processed_transcript"),"gene_id_short"]
  )
  
  #biotypes.vec = c("antisense","lincRNA","processed_transcript","all_biotypes","remove_antisense","lincProcTrans")
  #biotypes.vec = c("lincProcTrans", "lincRNA")
  
  
  
  ## EigenRank Setup
  sourceFile <- getFullPath("data/lncExprWithStats_transEachSample.tab")
  df <- readInTable(sourceFile)
  df <- df[which(df$averageExpr != 0),]   #df.1 = df.1[which(df.1$averageExpr != 0),]
  df[["gene_id_short"]] <-  sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
  doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  
  cols.list <- list(lpa=exprCols.lpa,lnpa=exprCols.lnpa,bothPullDowns=c(exprCols.lpa,exprCols.lnpa))
  
  lncFound.df <- getEnslist()
  lncFound.df <- unique(data.frame(ensembl_gene_id=lncFound.df$ensembl_gene_id,external_gene_id=lncFound.df$external_gene_id,gene_biotype=lncFound.df$gene_biotype))
  lncFound.df <- lncFound.df[which(lncFound.df$ensembl_gene_id %in% df[["gene_id_short"]]), ]
  
  df$label <- 0
  df[which(df$gene_id_short %in% lncFound.df$ensembl_gene_id),"label"] <- 1
  
    qc.df <- preProcessData()
  if (1 == 1){
    rows.list <- list()
    rows.list[["IDRlessthan0_01"]] <- qc.df[which(qc.df$IDR < 0.01),"gene_id_short"]
    rows.list[["IDRlessthan0_1"]]  <- qc.df[which(qc.df$IDR < 0.1),"gene_id_short"]
    rows.list[["IDRlessthan0_2"]]  <- qc.df[which(qc.df$IDR < 0.2),"gene_id_short"]
    rows.list[["IDRnotNA"]]        <- qc.df[which(!is.na(qc.df$IDR)),"gene_id_short"]
    rows.list[["allLncRNA"]]       <- qc.df[,"gene_id_short"]
    #lncGeneIdShorts.vec = c("IDRlessthan0_01", "IDRlessthan0_1","IDRlessthan0_2","allLncRNA")
    rows.list = lapply(rows.list, unique)
  }
  
  analysisPref = list(
    pca = TRUE,
    page = TRUE,
    logreg = TRUE,
    compare = TRUE)
#  analysisPref$compare = FALSE; analysisPref$pca = FALSE;analysisPref$page = FALSE; analysisPref$logreg = TRUE;
  
   print("starting run ... all libs are loaded...")
  
  #for(columns in c("lpa","lnpa","bothPullDowns")){
  #biotype <- "lincProcTrans";columns <- "lpa"; lncGroup <- "IDRlessthan0_1"
  for (biotype in biotypes.vec){
    lncRNA.biotype <- biotypes.list[[biotype]]
    basedir = paste(outdirFull, "/", biotype,"/", sep = "")
    for(columns in cols.vec){ 
      expr.cols <- cols.list[[columns]]
      for(lncGroup in lncGeneIdShorts.vec){
      
        expr.rows.group <- rows.list[[lncGroup]]
        expr.rows <- expr.rows.group[which(expr.rows.group %in% lncRNA.biotype )]
        
        expr.cols <- cols.list[[columns]]
        
        print(paste("starting",columns,lncGroup))
        outdir = paste(basedir,columns,"-",lncGroup,sep="")
        
        if(!file.exists(outdir)){dir.create(outdir)}
        outdir = paste(outdir,"/",sep="")
        print(outdir)
        
        local.df = df[which(df$gene_id_short %in% as.vector(expr.rows)),]
        local.df = local.df[which(!apply(local.df[expr.cols],1,function(x)sum(x)) == 0),]
        
        exportAsTable(df=local.df, file= paste(outdir,"comparison/",columns,"-",lncGroup,".tab",sep=""))
        if(sum(local.df$label) < 2){next}
        cat("sum(local.df$label",sum(local.df$label),lncGroup, columns,biotype)
        
        plotMsg = paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(local.df$label), "/", 
                        dim(local.df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ")      
        #    plotLncRNAComparisons <- function(lncDf=local.df,outdir=outdir,
        #                                      filename=paste(columns,biotype,sep="-"),cols=expr.cols ,titleMsg="",
        #                                      foundColword="",filebase=""){
        ## Log. Regression
        if (identical(analysisPref$logreg,TRUE)){
          print("logistic regression")
          
          bestCols.df <- data.frame(
            cols = expr.cols,
            score = sapply(expr.cols,function(x){mean(local.df[which(local.df$label == 1),x]) /mean(local.df[which(local.df$label == 0),x])}))
          topNCols = 10
          exprColsBestIndex = bestCols.df[order(bestCols.df$score,decreasing=TRUE)[1:topNCols],"cols"]
          exprLogReg = expr.cols[which(expr.cols %in% exprColsBestIndex)] 
          
          
          # logisticRegPcaExprData.R
          #runLogReg(lncDf=local.df,outdir=paste(outdir,"logReg",sep=""),
          #          cols=exprLogReg,iter=5,debug= FALSE,
          #          filebase=paste(columns,lncGroup,sep="-"),
          #          titleMsg=plotMsg)
          print("here")
          
          
          testGLMModel(outdir = paste(outdir,"logReg/",sep=""),
                       ratio=0.7,scaleDistro=FALSE,trials=10,msg=plotMsg,skipCustom=FALSE)
          
          #runLogRegIter(lncDf=local.df,cols=exprLogReg,outdir=paste(outdir,"logReg",sep=""),msg=plotMsg,trials=10)
          
          #runLogRegIter(lncDf=local.df,cols=exprLogReg,outdir=paste(outdir,"logReg",sep=""),msg=plotMsg,trials=1)
          
        } # end log reg
        if (identical(analysisPref$compare,TRUE)){
          print("comparison") 
          plotLncRNAComparisons(lncDf=local.df,outdir=paste(outdir,"/comparison/",sep=""),
                                filebase=paste(columns,lncGroup,sep="-"), cols=expr.cols,
                                titleMsg=plotMsg) #need label == 1
        }
        #lncDf=local.df;outdir=paste(outdir,"/comparison/",sep="");filebase=paste(columns,lncGroup,sep="-"); cols=expr.cols;titleMsg=plotMsg
        
        if (identical(analysisPref$page,TRUE)){
          print("page rank")
          # find fn in eigenRank_08272013.R
          #plotEigenVectorsDensity(lab=lab,nolab=nolab,outdir=outdir,
          #                        filename=paste(columns,biotype,sep="-"),cols=expr.cols ,titleMsg="")
          
          #plot w/o tissSpec == 1 ( what should by on axis is easier to cluster, w/o ts == 1 we can remove axis effects)
          #local.df <- df[which(df$gene_id_short %in% as.vector(expr.rows)),]
          
          nolab <- local.df[which(local.df$label == 0),]
          lab <- local.df[which(local.df$label == 1),]
          plotEigenVectorsDensity(lab=lab,nolab=nolab,outdir=paste(outdir,"/pageRank/",sep=""),
                                  filename=paste(columns,lncGroup,sep="-"),cols=expr.cols ,titleMsg=plotMsg)
          
          localts.df <- local.df[which(local.df$tissSpec != 1),]
          nolabts <- localts.df[which(localts.df$label == 0),]
          labts <- localts.df[which(localts.df$label == 1),]
          plotMsgTs = paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(localts.df$label), "/", 
                            dim(localts.df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ") 
          plotEigenVectorsDensity(lab=labts,nolab=nolabts,outdir=paste(outdir,"/pageRank/",sep=""),
                                  filename=paste("tissSpecNotOne",columns,lncGroup,sep="-"),cols=expr.cols ,titleMsg=plotMsgTs)
        } # end page
        
        ## PCA 
        if (identical(analysisPref$pca,TRUE)){
          print("PCA w/ outlier removal")
          
          pcaAnalysisRemoveOutliersSelectOutliers(lncDf=local.df,exprCols=which(colnames(local.df) %in% expr.cols),
                                                  foundColword="LncRNA",filebase=paste(columns,lncGroup,sep="-"),outDir=paste(outdir,"pca",sep=""),
                                                  titleMsg=plotMsg)
        }

      }
    }
    
    
  }
}


runLogRegTest <-function(outdir =  getFullPath("plots/fullAnalysisExperiment/")){
  
  bm73.file= getFullPath("data/lnc_biotype_ens73.tab")
  bm.df <- readInTable(file=bm73.file)
  bm.df$bm_biotype <- bm.df$gene_biotype
  bm.df$gene_id_short <- bm.df$ensembl_gene_id
  
  biotypes.list = list(antisense = bm.df[which(bm.df$bm_biotype == "antisense"),"gene_id_short"],
                       processed_transcript = bm.df[which(bm.df$bm_biotype == "processed_transcript"),"gene_id_short"],
                       lincRNA = bm.df[which(bm.df$bm_biotype == "lincRNA"),"gene_id_short"],
                       all_biotypes = bm.df[["gene_id_short"]],
                       remove_antisense =  bm.df[which(bm.df$bm_biotype != "antisense"),"gene_id_short"],
                       lincProcTrans = bm.df[which(bm.df$bm_biotype == "lincRNA" | bm.df$bm_biotype == "processed_transcript"),"gene_id_short"]
  )
  
  biotypes.vec = c("antisense","lincRNA","processed_transcript","all_biotypes","remove_antisense","lincProcTrans")
  
  
  ## EigenRank Setup
  sourceFile = getFullPath("data/lncExprWithStats_transEachSample.tab")
  df = readInTable(sourceFile)
  df = df[which(df$averageExpr != 0),]   #df.1 = df.1[which(df.1$averageExpr != 0),]
  df[["gene_id_short"]] =  sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
  doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  cols.list = list(lpa=exprCols.lpa,lnpa=exprCols.lnpa,bothPullDowns=c(exprCols.lpa,exprCols.lnpa))
  
  lncFound.df = getEnslist()
  lncFound.df = unique(data.frame(ensembl_gene_id=lncFound.df$ensembl_gene_id,external_gene_id=lncFound.df$external_gene_id,gene_biotype=lncFound.df$gene_biotype))
  lncFound.df = lncFound.df[which(lncFound.df$ensembl_gene_id %in% df[["gene_id_short"]]), ]
  
  df$label = 0
  df[which(df$gene_id_short %in% lncFound.df$ensembl_gene_id),"label"] = 1
  
  #biotype.df = readInTable(getFullPath("data/lnc_biotype.tab"))
  #df = merge(df,biotype.df[c("gene_id_short","bm_biotype")],by.x="gene_id_short",by.y="gene_id_short")  
  
  qc.df <- preProcessData()
  if (1 == 1){
    rows.list <- list()
    rows.list[["IDRlessthan0_01"]] <- qc.df[which(qc.df$IDR < 0.01),"gene_id_short"]
    rows.list[["IDRlessthan0_1"]]  <- qc.df[which(qc.df$IDR < 0.1),"gene_id_short"]
    rows.list[["IDRlessthan0_2"]]  <- qc.df[which(qc.df$IDR < 0.2),"gene_id_short"]
    rows.list[["IDRnotNA"]]        <- qc.df[which(!is.na(qc.df$IDR)),"gene_id_short"]
    rows.list[["allLncRNA"]]       <- qc.df[,"gene_id_short"]
    lncGeneIdShorts.vec = c("IDRlessthan0_01", "IDRlessthan0_1","IDRlessthan0_2","allLncRNA")
    rows.list = lapply(rows.list, unique)
  }
  
  analysisPref = list(
    pca = TRUE,
    page = TRUE,
    logreg = TRUE,
    compare = TRUE)
  analysisPref$pca = FALSE;analysisPref$page = FALSE; analysisPref$logreg = FALSE;
  
  print("starting run ... all libs are loaded...")
  
  biotype <- "lincProcTrans"
  lncRNA.biotype <- biotypes.list[[biotype]]
  basedir = paste(outdir, "/", biotype,"/", sep = "")
  columns = "lpa" 
  expr.cols <- cols.list[[columns]]
  lncGroup <- "IDRlessthan0_1"
  
  ## EigenRank
  expr.rows.group <- rows.list[[lncGroup]]
  expr.rows <- expr.rows.group[which(expr.rows.group %in% lncRNA.biotype )]
  
  expr.cols <- cols.list[[columns]]
  print(paste("starting",columns,lncGroup))
  outdir = paste(basedir,columns,"-",lncGroup,"/logRegressionTest",sep="")
  
  if(!file.exists(outdir)){dir.create(outdir)}
  outdir = paste(outdir,"/",sep="")
  
  local.df = df[which(df$gene_id_short %in% as.vector(expr.rows)),]
  local.df = local.df[which(!apply(local.df[expr.cols],1,function(x)sum(x)) == 0),]
  
  plotMsg = paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(local.df$label), "/", 
                  dim(local.df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ")      
  
  print("logistic regression")
  
  bestCols.df <- data.frame(
    cols = expr.cols,
    score = sapply(expr.cols,function(x){mean(local.df[which(local.df$label == 1),x]) /mean(local.df[which(local.df$label == 0),x])}))
  topNCols = length(expr.cols)
  exprColsBestIndex = bestCols.df[order(bestCols.df$score,decreasing=TRUE)[1:topNCols],"cols"]
  exprLogReg = expr.cols[which(expr.cols %in% exprColsBestIndex)] 
  
  
  #ratio value exper
  ratio.df <- expand.grid(mainEffects = c(TRUE),reg = c(FALSE),ratio=seq(0.5,0.1),run=1:10)
  ratio.df$predict <- NA
  for(i in 1:dim(ratio.df)[1]){
    print(i)
    
    line <- ratio.df[i,]
    output <- trainAndTestLogReg(lncDf=local.df,cols=exprLogReg, ratio = line$ratio, mainEffects=line$mainEffects,reg=line$reg)
    
    ratio.df$predict[i] <- output$R
    ratio.df$TPR[i] <- output$TPR
    ratio.df$TP[i] <- output$TP
    ratio.df$FP[i] <- output$FP
    ratio.df$FN[i] <- output$FN
    ratio.df$TN[i] <- output$TN
    ratio.df$AUC[i] <- output$AUC
  }
  exportAsTable(ratio.df, paste(outdir,"ratioExpr.tab",sep="/"))
  exportAsTable(local.df, paste(outdir,"fullExpr.tab",sep="/"))
  #lambda exper
  lambda.df <- as.data.frame(expand.grid(mainEffects = c(TRUE),ratio=0.7,reg = c(TRUE),lambda=c(0,10 ^ (seq(-2,8))),run=1:10))
  lambda.df$predict <- NA
  for(i in 1:dim(lambda.df)[1]){
    print(i)
    line <- lambda.df[i,]
    output <- trainAndTestLogReg(lncDf=local.df,cols=exprLogReg, ratio = line$ratio, mainEffects=line$mainEffects,reg=line$reg, lambda = line$lambda)
    lambda.df$predict[i] <- output$R
    lambda.df$TPR[i] <- output$TPR
    
    lambda.df$TP[i] <- output$TP
    lambda.df$FP[i] <- output$FP
    lambda.df$FN[i] <- output$FN
    lambda.df$TN[i] <- output$TN
    lambda.df$AUC[i] <- output$AUC
    
  }
  exportAsTable(lambda.df, paste(outdir,"lambdaExpr.tab",sep="/"))
  
}
#  runLogRegTest();analyzeLogRegTest()

analyzeLogRegTest <- function(outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/"),indir,msg){
  #print(outdir);print(indir);print(msg)
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  if(missing(indir)){
    indir = getFullPath("plots/fullAnalysisExperiment/lincProcTrans/lpa-IDRlessthan0_1/logRegressionTest/")
  }
  if (missing(msg)){
    
    #msg = paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(local.df$label), "/", 
    #                dim(local.df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ")      
    
    
    msg <- "lncRNA group = IDRlessthan0_1 ::: data cols = lpa \n( 42 / 1363 ) = (func/total lncRNA) ::: biotype = remove_antisense"
  }
  ratio.df <- readInTable(paste(indir,"ratioExpr.tab",sep="/"))
  ratio.df <- ratio.df[which(ratio.df$ratio < 1),]
  ratio.df <- within(ratio.df,{
    sens = TP / (TP + FN)
    prec = TP / (TP + FP)
    F1 = 2*(prec*sens)/(prec+sens)})
  ratddply.df <- ddply(ratio.df ,.(ratio), summarise, stdDev = sd(predict),mean=mean(predict), sensMean = mean(sens), sensSd = sd(sens), precMean=mean(prec), precSd = sd(prec))
  
  
  ggplot(ratio.df, aes(x=as.factor(ratio),y=AUC)) + geom_boxplot() +
    ylab("AUC = area under curve") +xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ ratio = train/total",msg, sep="\n"))
  ggsave(paste(outdir,"ratioAUC.pdf",sep=""))
  
  ggplot(ratio.df, aes(x=as.factor(ratio),y=sens)) + geom_boxplot() +
    ylab("Sensitivity") +xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ ratio = train/total",msg, sep="\n"))
  ggsave(paste(outdir,"ratioSensitivity.pdf",sep=""))
  
  ggplot(ratio.df, aes(x=as.factor(ratio),y=F1)) + geom_boxplot() +
    ylab("F1 = 2*(prec*sens)/(prec+sens) ") +xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ ratio = train/total",msg, sep="\n"))
  ggsave(paste(outdir,"ratioF1.pdf",sep=""))
  
  
  ggplot(ratio.df, aes(x=as.factor(ratio),y=prec)) + geom_boxplot()+
    ylab("Precision") +xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ ratio = train/total",msg, sep="\n"))
  ggsave(paste(outdir,"ratioPrecision.pdf",sep=""))
  
  ggplot(ratio.df, aes(x=as.factor(ratio),y=predict)) + geom_boxplot()+
    ylab("(TP + TN )/ total") +xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ ratio = train/total",msg, sep="\n"))
  ggsave(paste(outdir,"ratioExpr.pdf",sep=""))
  
  ggplot(ratio.df, aes(x=as.factor(ratio),y=TPR)) + geom_boxplot()+
    ylab("TP / P") + xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ ratio = train/total",msg, sep="\n"))
  ggsave(paste(outdir,"ratioExpr-TPoverP.pdf",sep=""))  
  
  sensPrec.df<-melt(ratio.df, id.var="ratio",measure.var=c("sens", "prec"))
  ggplot(sensPrec.df, aes(x=factor(ratio), y = value, fill=variable))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("ratio(test/total)") + ylab("paramater value from trials")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"ratioExpr-PrecSens.pdf",sep=""),height=7,width=9)
  
  ggplot(ratio.df,  aes(factor(ratio),predict))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("ratio(test/total)") + ylab("1 - test set error")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"ratioExpr-1-errorRate.pdf",sep=""),height=7,width=9)
  
  ggplot(ratio.df,  aes(factor(ratio),2*(prec*sens)/(prec+sens)))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("ratio(test/total)") + ylab("F1 = 2*(prec*sens)/(prec + sens) where range={0,1}")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n"))+ 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"ratioExpr-F1.pdf",sep=""),height=7,width=9)
  
  ggplot(ratio.df,  aes(factor(ratio),(TP + TN)/( TP + TN + FP + FN  )))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("ratio(test/total)") + ylab("accuracy = (TP + TN)/( TP + TN + FP + FN  )")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"ratioExpr-accuracy.pdf",sep=""),height=7,width=9)
  
  ggplot(ratio.df,  aes(factor(ratio),(TN)/(TN + FP)))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("ratio(test/total)") + ylab("true negative rate = (TN)/(TN + FP)")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"ratioExpr-TNR.pdf",sep=""),height=7,width=9)
  
  
  
  lambda.df <- readInTable(paste(indir,"lambdaExpr.tab",sep="/"))
  lambda.df <- within(lambda.df,{
    sens = TP / (TP + FN)
    prec = TP / (TP + FP)
    F1 =  2*(prec*sens)/(prec+sens)})
  
  lambdaddply.df <- ddply(lambda.df ,.(lambda), summarise,stdDev = sd(predict),mean=mean(predict), sensMean = mean(sens), sensSd = sd(sens), precMean=mean(prec), precSd = sd(prec))
  
  ggplot(lambda.df, aes(x=as.factor(lambda),y=AUC)) + geom_boxplot()+
    ylab("AUC = area under curve") +xlab("lambda value") +
    ggtitle(paste("logistic Regression: 10 trials w/  regularize lambda",msg, sep="\n"))
  ggsave(paste(outdir,"lambdaAUC.pdf",sep=""))
  
  
  ggplot(lambda.df, aes(x=as.factor(lambda),y=sens)) + geom_boxplot()+
    ylab("Sensitivity") +xlab("lambda value") +
    ggtitle(paste("logistic Regression: 10 trials w/  regularize lambda",msg, sep="\n"))
  ggsave(paste(outdir,"lambdaSensitivity.pdf",sep=""))
  
  ggplot(lambda.df, aes(x=as.factor(lambda),y=F1)) + geom_boxplot() +
    ylab("F1 = 2*(prec*sens)/(prec+sens) ") +xlab("ratio(test/total)") +
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n"))
  ggsave(paste(outdir,"lambdaF1.pdf",sep=""))
  
  
  ggplot(lambda.df, aes(x=as.factor(lambda),y=prec)) + geom_boxplot()+
    ylab("Precision") +xlab("lambda value") +
    ggtitle(paste("logistic Regression: 10 trials w/  regularize lambda",msg, sep="\n"))
  ggsave(paste(outdir,"lambdaPrecision.pdf",sep=""))
  
  ggplot(lambda.df, aes(x=as.factor(lambda),y=predict)) + geom_boxplot()+
    ylab("(TP + TN )/ total") + xlab("lambda value") +
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n"))
  ggsave(paste(outdir,"lambdaExpr.pdf",sep=""))
  
  ggplot(lambda.df, aes(x=as.factor(lambda),y=TPR)) + geom_boxplot()+
    ylab("TP / P") + xlab("lambda value") +
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n"))
  ggsave(paste(outdir,"lambdaExpr-TPR.pdf",sep=""))
  
  sensPrec.df<-melt(lambda.df, id.var="lambda",measure.var=c("sens", "prec"))
  ggplot(sensPrec.df, aes(x=factor(lambda), y = value, fill=variable))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("lambda value") + ylab("paramater value from trials")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"lambdaExpr-PrecSens.pdf",sep=""),height=7,width=9)
  
  #####
  ggplot(lambda.df,  aes(factor(lambda),predict))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("lambda value") + ylab("1 - test set error")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"lambdaExpr-1-errorRate.pdf",sep=""),height=7,width=9)
  
  ggplot(lambda.df,  aes(factor(lambda),2*(prec*sens)/(prec+sens)))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("lambda value") + ylab("F1 = 2*(prec*sens)/(prec + sens) where range={0,1}")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n"))+ 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"lambdaExpr-F1.pdf",sep=""),height=7,width=9)
  
  ggplot(lambda.df,  aes(factor(lambda),(TP + TN)/( TP + TN + FP + FN  )))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("lambda value") + ylab("accuracy = (TP + TN)/( TP + TN + FP + FN  )")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"lambdaExpr-accuracy.pdf",sep=""),height=7,width=9)
  
  ggplot(lambda.df,  aes(factor(lambda),(TN)/(TN + FP)))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("lambda value") + ylab("true negative rate = (TN)/(TN + FP)")+
    ggtitle(paste("logistic Regression: 10 trials w/ regularize lambda",msg, sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"lambdaExpr-TNR.pdf",sep=""),height=7,width=9)
  
}


testGLMforLogRegression <- function(df,ratio,scaleDistro=TRUE,skipCustom=FALSE){
  local.df <- df
  
  doubleCols = colnames(local.df)[as.vector(sapply(colnames(local.df),function(x)(typeof(local.df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  cols.list = list(lpa=exprCols.lpa,lnpa=exprCols.lnpa,bothPullDowns=c(exprCols.lpa,exprCols.lnpa))
  
  y <- local.df$label
  train <<- sample(seq_along(y), ratio * length(y))
  testFuncs <- sum(y[-train])
  y.train <- y[train]
  y.test  <- y[-train]
  y.freq0 <- length(which(y.train == 0))/length(y.train)
  y.freq1 <- length(which(y.train == 1))/length(y.train)
  weight.train <<- ifelse(y.train== 1, 1/y.freq1, 1/y.freq0)  
  ## scale training and test set seperately
  if( TRUE == scaleDistro ){
    local.df[train,exprCols.lpa] <- scale(local.df[train,exprCols.lpa], center=TRUE, scale = TRUE)
    local.df[-train,exprCols.lpa] <- scale(local.df[-train,exprCols.lpa], center=TRUE, scale = TRUE)
  }
  
  form.cross <- paste("label ~", do.call(paste, c(as.list(do.call(paste, c(expand.grid(exprCols.lpa, exprCols.lpa), sep=":"))), sep=" + ")))
  # use this formula
  print(paste(length(train),length(y), length(weight.train), dim(local.df[train,])[1]))
  glm.cross.weights <- glm(form.cross, data = local.df[train,], family=binomial, weights = weight.train) 
  glm.cross.weights.stats <- getStatsFromGlmModel(predict(glm.cross.weights,local.df[-train,], type="response"), y.test) 
  
  print("regular")
  glm.cross <- glm(form.cross, data = local.df, family=binomial,subset = train) 
  glm.cross.stats <- getStatsFromGlmModel(predict(glm.cross,local.df[-train,], type="response"), y.test) 
  print("     cross")  
  
  # log Reg with polynomial of degree two
  
  glm.poly2.weight <- glm(label ~ poly(HSMM.longPolyA, BJ.longPolyA, NHEK.longPolyA, HUVEC.longPolyA, MCF.7.longPolyA, K562.longPolyA, 
                                       AG04450.longPolyA, HMEC.longPolyA, HEPG2.longPolyA, SK.N.SH_RA.longPolyA, A549.longPolyA, H1.HESC.longPolyA, 
                                       NHLF.longPolyA, HELA.S3.longPolyA, GM12878.longPolyA, degree=2, raw=TRUE), data=local.df[train,], family = binomial, weights=weight.train)
  glm.poly2.weight.stats = getStatsFromGlmModel(predict(glm.poly2.weight,local.df[-train,], type="response"),      y.test)
  
  glm.poly2 <- glm(label ~ poly(HSMM.longPolyA, BJ.longPolyA, NHEK.longPolyA, HUVEC.longPolyA, MCF.7.longPolyA, K562.longPolyA, 
                                AG04450.longPolyA, HMEC.longPolyA, HEPG2.longPolyA, SK.N.SH_RA.longPolyA, A549.longPolyA, H1.HESC.longPolyA, 
                                NHLF.longPolyA, HELA.S3.longPolyA, GM12878.longPolyA, degree=2, raw=TRUE), data=local.df, family = binomial,subset=train)
  glm.poly2.stats = getStatsFromGlmModel(predict(glm.poly2,local.df[-train,], type="response"),      y.test)
  
  
  
  print("     poly")  
  glm.main.weight <- glm(label ~ HSMM.longPolyA + BJ.longPolyA + NHEK.longPolyA+ HUVEC.longPolyA+ MCF.7.longPolyA+ K562.longPolyA+ 
                           AG04450.longPolyA+ HMEC.longPolyA+ HEPG2.longPolyA+ SK.N.SH_RA.longPolyA+ A549.longPolyA+ H1.HESC.longPolyA+ 
                           NHLF.longPolyA+ HELA.S3.longPolyA+GM12878.longPolyA, data=local.df[train,], family = binomial, weights=weight.train)
  glm.main.weight.stats = getStatsFromGlmModel(predict(glm.main.weight,local.df[-train,], type="response"), y.test) 
  
  glm.main <- glm(label ~ HSMM.longPolyA + BJ.longPolyA + NHEK.longPolyA+ HUVEC.longPolyA+ MCF.7.longPolyA+ K562.longPolyA+ 
                    AG04450.longPolyA+ HMEC.longPolyA+ HEPG2.longPolyA+ SK.N.SH_RA.longPolyA+ A549.longPolyA+ H1.HESC.longPolyA+ 
                    NHLF.longPolyA+ HELA.S3.longPolyA+GM12878.longPolyA, data=local.df, family = binomial,subset=train)
  glm.main.stats = getStatsFromGlmModel(predict(glm.main,local.df[-train,], type="response"), y.test) 
  print("     main")  
  
  
  
  form.circ <-paste("label ~ ",paste("I(", do.call(paste,c(as.list(paste0(exprCols.lpa,paste0("*",exprCols.lpa))),sep=") + I("  )), ")"))
  
  glm.circ.weights <- glm(form.circ, data = local.df[train,], family=binomial,weights = weight.train) 
  glm.circ.weights.stats <- getStatsFromGlmModel(predict(glm.circ.weights,local.df[-train,], type="response"), y.test) 
  
  
  glm.circ <- glm(form.circ, data = local.df, family=binomial,subset = train) 
  glm.circ.stats <- getStatsFromGlmModel(predict(glm.circ,local.df[-train,], type="response"), y.test) 
  print("     circ")  
  
  
  lda.main <- lda(label ~ HSMM.longPolyA + BJ.longPolyA + NHEK.longPolyA+ HUVEC.longPolyA+ MCF.7.longPolyA+ K562.longPolyA+ 
                    AG04450.longPolyA+ HMEC.longPolyA+ HEPG2.longPolyA+ SK.N.SH_RA.longPolyA+ A549.longPolyA+ H1.HESC.longPolyA+ 
                    NHLF.longPolyA+ HELA.S3.longPolyA+GM12878.longPolyA, data=local.df, subset=train)
  lda.main.stats = getStatsFromGlmModel(predict(lda.main,local.df[-train,])$x, y.test) 
  print("     lda")  
  
  qda.main <- qda(label ~ HSMM.longPolyA + BJ.longPolyA + NHEK.longPolyA+ HUVEC.longPolyA+ MCF.7.longPolyA+ K562.longPolyA+ 
                    AG04450.longPolyA+ HMEC.longPolyA+ HEPG2.longPolyA+ SK.N.SH_RA.longPolyA+ A549.longPolyA+ H1.HESC.longPolyA+ 
                    NHLF.longPolyA+ HELA.S3.longPolyA+GM12878.longPolyA, data=local.df, subset=train)
  qda.main.stats = getStatsFromGlmModel(predict(qda.main,local.df[-train,], type="response")$posterior[,2], y.test) 
  print("     qda")  
  
  #KNN
  train.X <- as.matrix(local.df[train,exprCols.lpa])
  test.X <- as.matrix(local.df[-train,exprCols.lpa])
  direction <- local.df$label
  train.Direction <- direction[train]
  
  knn.1.stats <- getStatsFromGlmModel(knn(train.X, test.X, train.Direction, k = 1), y.test, knn=TRUE)
  
  knn.3.stats <- getStatsFromGlmModel( knn(train.X, test.X, train.Direction, k = 3), y.test, knn=TRUE)
  
  
  knn.5.stats <- getStatsFromGlmModel(knn(train.X, test.X, train.Direction, k = 5), y.test, knn=TRUE)
  print("     knn")  
  
  grid = 10 ^ seq(10, -2, length = 100)
  ridge.weight.stats <- getStatsFromGlmModel( predict(glmnet(train.X,y.train,family="binomial", alpha=0, lambda=grid, weights = weight.train,standardize=scaleDistro),
                                                      s=cv.glmnet(train.X,y.train,family="binomial",type.measure="auc", alpha=0, weights=weight.train,standardize=scaleDistro)$lambda.min,
                                                      newx=test.X, type="response")  , 
                                              y.test)
  print("     ridge.weight.stats") 
  
  lasso.weight.stats <- getStatsFromGlmModel(  predict( glmnet(train.X,y.train,family="binomial", alpha=1, lambda=grid, weights = weight.train,standardize=scaleDistro),
                                                        s=cv.glmnet(train.X,y.train,family="binomial",type.measure="auc", alpha=1, weights=weight.train,standardize=scaleDistro)$lambda.min,
                                                        newx=test.X, type="response"), 
                                               y.test)
  print("     lasso.weight.stats") 
  lasso.stats <- getStatsFromGlmModel( predict(glmnet(train.X,y.train,family="binomial", alpha=1, lambda=grid, weights = rep(1,length(y.train)),standardize=scaleDistro),
                                               s=cv.glmnet(train.X,y.train,family="binomial",type.measure="auc", alpha=1, weights=rep(1,length(y.train)),standardize=scaleDistro)$lambda.min,
                                               newx=test.X, type="response"),
                                                y.test)
  print("     lasso.stats") 
  
  ridge.stats <- getStatsFromGlmModel( predict( object = glmnet(train.X,y.train,family="binomial", alpha=0, lambda=grid, weights = rep(1,length(y.train)),standardize=scaleDistro),
                                                s = cv.glmnet(train.X,y.train,family="binomial",type.measure="auc", alpha=0, weights=rep(1,length(y.train)),standardize=scaleDistro)$lambda.min,
                                                newx = test.X, type="response"), 
                                        y.test)
  print("     ridge/ridge.stats")  
  
  
  
  ellipse.stats <- trainAndTestLogReg(lncDf=local.df,cols=exprCols.lpa, ratio = ratio,
                                      mainEffects=TRUE,reg=TRUE, lambda = 0,glmTypeOutput=TRUE )
  print("      ellipse algo stats")
  
  
  df.out <-  rbind(glm.poly2.weight.stats,glm.poly2.stats,glm.main.weight.stats,glm.main.stats,
                   glm.cross.weights.stats,glm.cross.stats,lda.main.stats,qda.main.stats,
                   knn.1.stats,knn.3.stats,knn.5.stats,glm.circ.weights.stats,glm.circ.stats,
                   ridge.weight.stats,ridge.stats,lasso.weight.stats,lasso.stats,ellipse.stats)
  
  df.out$expr <- c("glm.poly2.weight.stats","glm.poly2.stats","glm.main.weight.stats","glm.main.stats",
                   "glm.cross.weights.stats","glm.cross.stats","lda.main.stats","qda.main.stats",
                   "knn.1.stats","knn.3.stats","knn.5.stats","glm.circ.weights.stats","glm.circ.stats",
                   "ridge.weight.stats","ridge.stats","lasso.weight.stats","lasso.stats","ellipse.stats")
  rownames(df.out) <- NULL
  df.out$testFunctionalLncs <- testFuncs 
  df.out$testSetSize <- length(train)
  df.out
  
}


# rm(list=unlist(ls())); source("~/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/fullAnalysis-Aug2013.R"); testGLMModelForRatios()

testGLMModelForRatios <- function(){
  odir <- getFullPath("plots/fullAnalysisExperiment/test/logReg/glmTest-wEllipseAlgo")
  
  ratio1.odr <- paste(odir,"/", sep="")
  testGLMModel(outdir=ratio1.odr, ratio=0.7)
  
  ratio1.odr <- paste(odir,"scale/", sep="")
  testGLMModel(outdir=ratio1.odr, ratio=0.7,scaleDistro=TRUE)
}

# rm(list=unlist(ls())); source('~/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/fullAnalysis-Aug2013.R"); testGLMModel()
testGLMModel <- function(outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/glmTest-wEllipseAlgo/"),
                         dir,ratio=0.7,scaleDistro=FALSE,trials=30,msg,df,skipCustom=FALSE){
  if (!file.exists(outdir)){
    dir.create(outdir)}
 

  if(missing(df)){
    df <- readInTable(getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/fullExpr.tab"))  
  }
  
  if (missing(msg)){
    biotype <- "lincProcTrans";colmuns<-"lpa";lncGroup <- "IDRlessthan0_1";
    plotMsg = paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(df$label), "/", 
                    dim(df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ")      
  }
  
  
  exprFile <- paste(outdir,"glmExpr",trials,"-",ratio,".tab",sep="")
  test = TRUE
  if (!file.exists(exprFile) || test == TRUE){
    print(paste("ratio = ", ratio, "trial =  1"))
    df.accum <- testGLMforLogRegression(df,ratio,scaleDistro)
    df.accum$trial <- 1
    for(i in 2:trials){
      print(paste("trial = " ,i,sep=""))
     
      tryCatch({df.out <- testGLMforLogRegression(df,ratio,scaleDistro,skipCustom=skipCustom)
      df.out$trial <- i
      df.accum <- rbind(df.out, df.accum)},
      error = function(e) print("error in Log Reg routine"), finally=print(""))
    }
    exportAsTable(df.accum, exprFile)
    local.df <- df.accum
  } else{
    local.df <- readInTable(exprFile)  
  }
  
  summary.df <- ddply(local.df, .(expr), summarize, precMean= mean(prec,na.rm=T),sensMean= mean(sens,na.rm=T))
  sensPrec.df<-melt(local.df, id.var="expr",measure.var=c("sens", "prec"))
  biotype <- "lincProcTrans";colmuns<-"lpa";lncGroup <- "IDRlessthan0_1";
  msg <- paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(df$label), "/", 
                  dim(df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ") 
  msg1 <- paste("traning set Size = ", local.df$testSetSize[1], ",trials=", max(local.df$trial),"ratio = ", ratio, "scale = ", scaleDistro,sep=" " )
  
  ggplot(sensPrec.df, aes(x=factor(expr), y = value, fill=variable))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("learning algorithm") + ylab("paramater value from trials")+
    ggtitle(paste("test set error for classification algorithm",msg,msg1,sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"glmModelComparison-PrecSens.pdf",sep=""),height=7,width=9)
  
  ggplot(local.df,  aes(expr,1-errorRate))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("learning algorithm") + ylab("1 - test set error")+
    ggtitle(paste("test set error for classification algorithm",msg,msg1,sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"glmModelComparison-1-errorRate.pdf",sep=""),height=7,width=9)
  
  ggplot(local.df,  aes(expr,AUC))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("learning algorithm") + ylab("AUC")+
    ggtitle(paste("test set error for classification algorithm",msg,msg1,sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"glmModelComparison-AUC.pdf",sep=""),height=7,width=9)
  
  
  ggplot(local.df,  aes(expr,2*(prec*sens)/(prec+sens)))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("learning algorithm") + ylab("F1 = 2*(prec*sens)/(prec + sens) where range={0,1}")+
    ggtitle(paste("test set error for classification algorithm",msg,msg1,sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"glmModelComparison-F1.pdf",sep=""),height=7,width=9)
  
  ggplot(local.df,  aes(expr,(TP + TN)/( TP + TN + FP + FN  )))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("learning algorithm") + ylab("accuracy = (TP + TN)/( TP + TN + FP + FN  )")+
    ggtitle(paste("test set error for classification algorithm",msg,msg1,sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"glmModelComparison-accuracy.pdf",sep=""),height=7,width=9)
  
  ggplot(local.df,  aes(expr,(TN)/(TN + FP)))+geom_boxplot() + 
    coord_flip() + theme_bw()+ xlab("learning algorithm") + ylab("true negative rate = (TN)/(TN + FP)")+
    ggtitle(paste("test set error for classification algorithm",msg,msg1,sep="\n")) + 
    theme(panel.grid.major.x= element_line(colour = "grey"))+
    theme(panel.grid.major.y= element_line(colour = "black")) 
  ggsave(file=paste(outdir,"glmModelComparison-TNR.pdf",sep=""),height=7,width=9)
  
}



getStatsFromGlmModel <- function(probs, y,knn=FALSE){
  
  if (TRUE == knn){ 
    pred <- as.numeric(probs) - 1 
  } else {
    pred <- rep(0,length(probs))
    pred[which(probs > 0.5)] <- 1
  }
  
  correct <- (pred == y)
  poly2 <- data.frame(trial=-1)
  poly2$TP <- length(which(correct & y ==1))
  poly2$TN <- length(which(correct & y ==0))  
  poly2$FP <- length(which(!correct & y ==0))  
  poly2$FN <- length(which(!correct & y ==1))  
  poly2$prec <- with(poly2, TP / (TP + FP))
  poly2$sens <- with(poly2, TP / (TP + FN))
  poly2$errorRate <-  1 - sum(correct)/length(correct)
  if (TRUE == knn){ 
    poly2$AUC <- 0
  } else {
    poly2$AUC <- calcAUC(prob=probs, label=y)
  }
  poly2
}


plotTissSpecVaverageExpr <- function(outdir = getFullPath("plots/lncCompare/globalAnalysisExpr/")){
  if(!file.exists(outdir)){dir.create(outdir)}
  
  sourceFile = getFullPath("data/lncExprWithStats_transEachSample.tab")
  df = readInTable(sourceFile)
  df = df[which(df$averageExpr != 0),]   #df.1 = df.1[which(df.1$averageExpr != 0),]
  df[["gene_id_short"]] =  sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
  doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  
  lncFound.df = getEnslist()
  lncFound.df = unique(data.frame(ensembl_gene_id=lncFound.df$ensembl_gene_id,external_gene_id=lncFound.df$external_gene_id,gene_biotype=lncFound.df$gene_biotype))
  lncFound.df = lncFound.df[which(lncFound.df$ensembl_gene_id %in% df[["gene_id_short"]]), ]
  
  df$label = 0
  df[which(df$gene_id_short %in% lncFound.df$ensembl_gene_id),"label"] = 1
  
  biotype.df = readInTable(getFullPath("data/lnc_biotype.tab"))
  df = merge(df,biotype.df[c("gene_id_short","bm_biotype")],by.x="gene_id_short",by.y="gene_id_short")
  
  
  df.bigThree <- df[which(df$bm_biotype %in% c("antisense", "lincRNA", "processed_transcript")),]
  
  g1  <- ggplot(df.bigThree, aes(,x=averageExpr, y = tissSpec, color=factor(bm_biotype))) + 
    geom_point() + 
    theme_bw() +
    ggtitle("average Expr vs. tissue specificity\nlncRNA")
  g1
  ggsave(file=paste(outdir,"aveExprVtissSpec.pdf",sep="/"),height=6,width=8)
  
  g1.size <- ggplot(df.bigThree, aes(,x=averageExpr, y = tissSpec, color=factor(bm_biotype), size = label + 1)) + 
    geom_point() + 
    theme_bw() +
    ggtitle("average Expr vs. tissue specificity\nlncRNA\nfunctional examples are larger")+
    scale_size(range=c(1,3)) 
  g1.size
  ggsave(file=paste(outdir,"aveExprVtissSpec_func=Size.pdf",sep="/"),height=6,width=8)
  
  g1 + xlim(c(0,quantile(df$averageExpr, 0.95)))
  ggsave(file=paste(outdir,"aveExprVtissSpec_zoom.pdf",sep="/"),height=6,width=8)
  
  g1.size + xlim(c(0,quantile(df$averageExpr, 0.95)))
  ggsave(file=paste(outdir,"aveExprVtissSpec_func=Size_zoom.pdf",sep="/"),height=6,width=8)
  
  
  
  df$tissSpecUnity = ifelse(df$tissSpec == 1, 1, 0)
  N.tissSpecUnity = sum(df$tissSpecUnity)
  N.tissSpecNotUnity = dim(df)[1] - N.tissSpecUnity
  
  ggplot(df[which(df$tissSpecUnity == 1),], aes(x=averageExpr))+geom_density()+
    theme_bw()+
    ggtitle(paste("lncRNA where tissSpec == 1","\nN=",N.tissSpecUnity ,sep=""))
  ggsave(file=paste(outdir,"tissSpec=1_RNAexprDistro.pdf",sep="/"),height=6,width=8)
  
  ggplot(df, aes(x = averageExpr, fill = factor(tissSpecUnity)))+geom_density(alpha=I(0.4)) +
    theme_bw()+
    xlim(0,.003) +
    ggtitle(paste("lncRNA","\ntissSpec == 1 : ",N.tissSpecUnity, "\ntissSpec != 1 : ",N.tissSpecNotUnity,sep=""))
  ggsave(file=paste(outdir,"tissSpec=1_RNAexprDistro_vsRest.pdf",sep="/"),height=6,width=8)
  
  
}



#runLogRegPac


testRegByPackage <- function(outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/LiblineaR/")){
  trials <- 20
  
  if (!file.exists(outdir)){
    dir.create(outdir)}
  msg <- "lncRNA group = IDRlessthan0_1 ::: data cols = lpa \n( 46 / 1954 ) = (func/total lncRNA) ::: biotype = lincProcTrans"
  msg <- paste(msg, "\n trial = ", trials)
  
  df <- readInTable(getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/fullExpr.tab"))
  
  exprCols.lpa  = getLpaColnames(df)
  
  x=df[,exprCols.lpa]
  y=factor(df$label)
  ratio <- 0.7
  
  t = 0 
  tryTypes = c(0,6,7)
  tryCosts = rev(c(10^seq(-5,10)))
  
  
  for(trial in 1:trials){
    train = sample(seq_along(y),length(y) * ratio )
    xTrain=x[train,]
    s = scale(xTrain,center=TRUE,scale=TRUE)
    xTest=x[-train,]
    yTrain=y[train]
    yTrain.freq1 = sum(yTrain == "1")/length(yTrain)
    yTrain.freq0 = 1 - yTrain.freq1
    yTest = y[-train]
    
    
    weight <- TRUE;ttype <-3
    for (co in tryCosts){
        
      if (TRUE == weight){
        acc = LiblineaR(data=s, labels=yTrain, type = ttype, cost = co, bias=TRUE,  verbose=FALSE, wi=list("0"=1/yTrain.freq0, "1"=1/yTrain.freq0))
      } else {
        acc = LiblineaR(data=s, labels=yTrain, type = ttype, cost = co, bias=TRUE,  verbose=FALSE)
      }
      
      cat("Results for C=", co, " : ", trial, " type \n", sep="")
      
      inLineTest <- TRUE
      if(TRUE == inLineTest){
        s2 = scale(xTest, attr(s, "scaled:center"), attr(s, "scaled:scale"))
        pr = FALSE
        if (ttype == 0 || ttype == 7) pr =TRUE
        p=predict(acc,s2,proba=pr, decisionValues=TRUE)
        res = table(p$predictions, yTest)
        if(dim(res)[1] > 1) print(res)
        df.out <- getStatsFromGlmModel(as.numeric(p$predictions) , yTest)
        df.out$cost <- co
        df.out$trial <- trial
        if (trial == 1 && co == max(tryCosts) ){
          df.accum <- df.out
        } else {
          df.accum <- rbind(df.accum, df.out)
          
        }
      }
    }
  }
  ggplot(melt(df.accum, id.var = "cost", measure.var = c("sens", "prec")), aes(x=factor(cost),y=value,fill=variable))+ 
    geom_boxplot() + theme_bw() + 
    ggtitle(paste("L2 regularized L1-loss support vector classification\n",msg,sep=""))
  ggsave(paste(outdir, "L2-regularized-L1-loss-SVC_outputOverC.pdf", sep=""))
  
  ggplot(melt(df.accum, id.var = "cost", measure.var = c("AUC")), aes(x=factor(cost),y=value,fill=variable))+ 
    geom_boxplot() + theme_bw() + 
    ggtitle(paste("L2 regularized L1-loss support vector classification\n",msg,sep=""))
  ggsave(paste(outdir, "L2-regularized-L1-loss-SVC_outputOverC-AUC.pdf", sep=""))
  
  ggplot(melt(df.accum, id.var = "cost", measure.var = c("errorRate")), aes(x=factor(cost),y=value,fill=variable))+ 
    geom_boxplot() + theme_bw() + 
    ggtitle(paste("L2 regularized L1-loss support vector classification\n",msg,sep=""))
  ggsave(paste(outdir, "L2-regularized-L1-loss-SVC_outputOverC-errorRate.pdf", sep=""))
    
  exportAsTable(df=df.accum, file = paste(outdir, "L2-regularized-L1-loss-SVC.tab", sep=""))
  #BCR = mean(c(res[1,1]/sum(res[,1]), res[2,2]/sum(res[,2]), res[3,3]/sum(res[,3])  ))
  
}

# calcAUC <- function(prob, label){
#   AUC <-NA
#   tryCatch({
#   pre = prediction(predictions=prob, labels=label)
#   per = performance(pre, "tpr", "fpr")
#   AUC = (performance(pre, "auc"))@y.values[[1]] 
#   })
#   
#   AUC
# }

createROCcurve <- function(prob, label,outfile,title=""){
  pre = prediction(predictions=prob, labels=label)
  per = performance(pre, "tpr", "fpr")
  AUC = (performance(pre, "auc"))@y.values[[1]] 
  subtitle=paste("AUC:", AUC)
  pdf(outfile)
  plot(per,main=paste("ROC Curve",title,sep="\n"), xlab="True Positive Rate",ylab="False Positive Rate")
  dev.off()
}


testRegByGLMNET <- function(outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/glmnetCV/"), weight=TRUE){
  trials <- 10
  
  if (!file.exists(outdir)){
    dir.create(outdir)}
  msg <- "lncRNA group = IDRlessthan0_1 ::: data cols = lpa \n( 46 / 1954 ) = (func/total lncRNA) ::: biotype = lincProcTrans"
  msg <- paste(msg, "\n trial = ", trials)
  
  df <- readInTable(getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/fullExpr.tab"))
  doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  
  x=df[,exprCols.lpa]
  y=factor(df$label)
  ratio <- 0.7
  
  for(trial in 1:trials){
    
    train = sample(seq_along(y),length(y) * ratio )
    xTrain=x[train,]
    xTest=x[-train,]
    yTrain=y[train]
    yTrain.freq1 = sum(yTrain == "1")/length(yTrain)
    yTrain.freq0 = 1 - yTrain.freq1
    yTest = y[-train]
    if(TRUE == weight){
      weight.train <- ifelse(yTrain == 1, 1/yTrain.freq1, 1/yTrain.freq0)
    } else {
      weight.train <- rep(1,length(yTrain))
    }
    grid = 10 ^ seq(10, -2, length = 100)
    
    ridge.mod <-   glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial", alpha=0, lambda=grid, weights = rep(1,length(yTrain)))
    ridge <- cv.glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial",type.measure="auc", alpha=0, weights=rep(1,length(yTrain)))
    ridge.predict <- predict(ridge.mod,s=ridge$lambda.min,newx=as.matrix(xTest), type="response")  
    ridge.stats <- getStatsFromGlmModel( ridge.predict, yTest)
    ridge.stats$type <- "ridge"
    
    lasso.mod <-   glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial", alpha=1, lambda=grid, weights = rep(1,length(yTrain)))
    lasso <- cv.glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial",type.measure="auc", alpha=1, weights=rep(1,length(yTrain)))
    lasso.predict <- predict(lasso.mod,s=lasso$lambda.min,newx=as.matrix(xTest), type="response")  
    lasso.stats <- getStatsFromGlmModel( lasso.predict, yTest)
    lasso.stats$type <- "lasso"
    
    ridge.weight.mod <-   glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial", alpha=0, lambda=grid, weights = weight.train)
    ridge.weight <- cv.glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial",type.measure="auc", alpha=0, weights=weight.train)
    ridge.weight.predict <- predict(ridge.weight.mod,s=ridge.weight$lambda.min,newx=as.matrix(xTest), type="response")  
    ridge.weight.stats <- getStatsFromGlmModel( ridge.weight.predict, yTest)
    ridge.weight.stats$type <- "ridge.weight"
    
    lasso.weight.mod <-   glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial", alpha=1, lambda=grid, weights = weight.train)
    lasso.weight <- cv.glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial",type.measure="auc", alpha=1, weights=weight.train)
    lasso.weight.predict <- predict(lasso.weight.mod,s=lasso.weight$lambda.min,newx=as.matrix(xTest), type="response")  
    lasso.weight.stats <- getStatsFromGlmModel( lasso.weight.predict, yTest)
    lasso.weight.stats$type <- "lasso.weight"
    
    
    df.out <- rbind(lasso.stats,lasso.weight.stats,ridge.stats,ridge.weight.stats)
    df.out$trial = trial
    if (trial == 1){
      df.accum <- df.out
    } else {
      df.accum <- rbind(df.out, df.accum)
    }
    
  } # end iteration over trials
  ggplot(melt(df.accum, id.var = "type", measure.var=c("prec", "sens", "errorRate")), aes(type, value, fill=variable))  + 
    geom_boxplot() +
    ggtitle(paste("ridge and lasso regularized linear Regression\nLambda calculated by CV, final values from indepenet set(0.3 total)", msg, sep="\n"))
}    

plotAUC <- function(prob, label){
  pre = prediction(predictions=prob, labels=label)
  per = performance(pre, "tpr", "fpr")
  plot(perf,colorize=TRUE)
  
}

processPredicitonForROC <- function(df, label,exprCols){
  df$dist <- df[[label]]
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
    label = withinSubset
  })
  dfo.stats[c("predict","totalMembers","P","N","TP","FP","TN","FN","sens","prec","TPR","FPR","FNR","TNR","Fscore","accuracy","errorRate","dist","label")]
}

#Jan 8, 2013
#rm(list=unlist(ls()));source('~/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/fullAnalysis-Aug2013.R');runRidgeRegression()
runRidgeRegression <- function(outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/"), weight=TRUE){
  trials <- 200
  
  if (!file.exists(outdir)){dir.create(outdir)}
  
    
  msg <- "lncRNA group = IDRlessthan0_1 ::: data cols = lpa \n ::: biotype = lincProcTrans"
  msg <- paste(msg, "\n trial = ", trials)
  
  df <- read.table(getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/fullExpr.tab"),stringsAsFactors=FALSE,header=TRUE)
  doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  gene_id_predict.df <- as.data.frame(df[c("gene_id")])
  gene_id_predict.df$count <- 0
  x=df[,exprCols.lpa]
  y=factor(df$label)
  ratio <- 0.7
  
  for(trial in 1:trials){
    train = sample(seq_along(y),length(y) * ratio )
    xTrain=x[train,]
    xTest=x[-train,]
    xTest.index = seq_along(df$gene_id_short)[-1 * train]
    xTest.lncRnaName <- df$gene_id[-train]
    yTrain=y[train]
    yTrain.freq1 = sum(yTrain == "1")/length(yTrain)
    yTrain.freq0 = 1 - yTrain.freq1
    yTest = y[-train]
    
      weight.train <- ifelse(yTrain == 1, 1/yTrain.freq1, 1/yTrain.freq0)

    grid = 10 ^ seq(10, -2, length = 100)
    
       
    ridge.weight.mod <-   glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial", alpha=0, lambda=grid, weights = weight.train)
    ridge.weight <- cv.glmnet(as.matrix(xTrain),as.numeric(yTrain)-1,family="binomial",type.measure="auc", alpha=0, weights=weight.train)
    ridge.weight.predict <- predict(ridge.weight.mod,s=ridge.weight$lambda.min,newx=as.matrix(xTest), type="response")  
    ridge.weight.stats <- getStatsFromGlmModel( ridge.weight.predict, yTest)
    ridge.weight.stats$type <- "ridge.weight"
    print(paste(trial, "-----"))
    
    df.test <- data.frame(ridge.weight.predict)
    colnames(df.test) <- "probs"
    df.test$withinSubset <- as.numeric(yTest) - 1
    df.tpr <- processPredicitonForROC(df.test, label="probs") # "probs" -> "dist" in df colnames
   # print(paste(trial, "-----"))
    
    #ggplot(df.tpr, aes(FPR,TPR))+geom_abline(slope=1)+geom_smooth()+theme_bw()+xlim(0,1)+ylim(0,1)+geom_point()
    
    df.tpr$trial <- trial
    df.tpr$probs <- df.tpr$dist
    
    #print(dim(df.tpr))
    #print(length(xTest.lncRnaName))
    df.tpr$lncRnaName <- xTest.lncRnaName
    ridge.weight.stats$trial <- trial
    if (trial == 1){
      df.accum <- ridge.weight.stats
      df.accum.tpr <- df.tpr
    } else {
      df.accum <- rbind(df.accum, ridge.weight.stats)
      df.accum.tpr <- rbind(df.accum.tpr, df.tpr)
    }
  #  print(paste(trial, "-----"))
    
    #count the predicted correct lncRNA:
    gene_id_predict.df$count = gene_id_predict.df$count + ifelse(gene_id_predict.df$gene_id %in% xTest.lncRnaName[which(df.test$probs > 0.5)],1,0)
    
  } # end iteration over trials
  
  exportAsTable(df=df.accum, file=paste(outdir,"ridgeWeightStats-summary.tab",sep=""))
  exportAsTable(df=df.accum.tpr, file=paste(outdir,"ridgeWeightStats-ROC.tab",sep=""))
  exportAsTable(df=gene_id_predict.df, file=paste(outdir,"genePredictionCounts.tab",sep=""))
  
  ggplot(melt(df.accum, id.var = "type", measure.var=c("prec", "sens", "errorRate")), aes(type, value, fill=variable))  + 
    geom_boxplot() +
    ggtitle(paste("ridge and lasso regularized linear Regression\nLambda calculated by CV, final values from indepenet set(0.3 total)", msg, sep="\n"))
}    


#Jan 9, 2013 
analyzeRidgeReg <- function(outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/") ){
  sourceFile = getFullPath("data/lncExprWithStats_transEachSample.tab")
  lncGene.df = readInTable(sourceFile)
  
  
  ridgeStats.df <- readInTable(file=paste(outdir,"ridgeWeightStats-summary.tab",sep=""))
  ridgeTPR.df <- readInTable(file=paste(outdir,"ridgeWeightStats-ROC.tab",sep=""))
# ridgeTPR.df "probs" ~~ "dist" in df colnames
  ridgeGenePred.df <- readInTable(file=paste(outdir,"genePredictionCounts.tab",sep=""))
  lncFound.df = getEnslist()
  
  #find a list of lncRNA identified by ridge regression
  lncRidge.gene_id <- as.vector(sapply(ridgeGenePred.df[which(ridgeGenePred.df$count > 0), "gene_id"],
                                   function(x)unlist(strsplit(x, "\\."))[1]))
  recoveredFunc.df <- unique(lncFound.df[which(lncFound.df$ensembl_gene_id %in% lncRidge.gene_id),c("ensembl_gene_id","external_gene_id")])
  allFound.df <- getExternalGeneIdFromEnsGene(lncRidge.gene_id,attr=c("ensembl_gene_id","hgnc_symbol","gene_biotype"))
  allFound.df$handle <- ifelse(allFound.df$hgnc_symbol == "",allFound.df$ensembl_gene_id,allFound.df$hgnc_symbol)
  exportAsTable(df=allFound.df, file=paste(outdir,"foundGeneByRidgeRegression.tab",sep=""))
 
  ## explore biotypes further
  # "protein coding" biotype? 
  lncRidge.pc <- allFound.df[which(allFound.df$gene_biotype == "protein_coding"),"ensembl_gene_id"]
  lncGene.df[which(lncGene.df$lncRnaName %in% lncRidge.pc),]
  
  ## analyze the stats from each of the 200 trials
  # dim(ridgeTPR.df) => "[1] 117400     20"
  ggplot(ridgeStats.df, aes(AUC))+geom_density(alpha=I(0.4))+theme_bw()+xlim(0,1) + 
    ggtitle("(nonAS, IDR<0.1)Ridge Regression(weighted); AUC")
  ggsave(paste(outdir,"ridgeReg_AUC.pdf",sep=""))
  
  ggplot(ridgeStats.df, aes(errorRate))+geom_density(alpha=I(0.4))+theme_bw()+xlim(0,1) + 
    ggtitle("(nonAS, IDR<0.1)Ridge Regression(weighted); AUC")
  ggsave(paste(outdir,"ridgeReg_errorRate.pdf",sep=""))
  
  ## get the best(AUC) trial...
  # for median trial, we need an odd number, let's just add one(max AUC) 

  # plot ROC curves for median and max trial...
  best.trial <- ridgeStats.df[which(ridgeStats.df$AUC == max(ridgeStats.df$AUC)),"trial"]
  df.tpr.best <- ridgeTPR.df[which(ridgeTPR.df$trial == best.trial),]
  ggplot(df.tpr.best, aes(FPR,TPR))+geom_abline(slope=1)+ geom_smooth() + 
    theme_bw()+xlim(0,1)+ylim(0,1) + geom_point() +
    ggtitle("(nonAS, IDR<0.1)Ridge Reg(weighted) : best trial PR curve")
  ggsave(paste(outdir,"ridgeReg_ROC_best.pdf",sep=""))
  
  median.trial <- ridgeStats.df[which(ridgeStats.df$AUC == median(c(ridgeStats.df$AUC,1))),"trial"]
  df.tpr.median <- ridgeTPR.df[which(ridgeTPR.df$trial == median.trial),]
  ggplot(df.tpr.median, aes(FPR,TPR))+geom_abline(slope=1) + geom_smooth() + 
    theme_bw()+xlim(0,1)+ylim(0,1)+geom_point()+
    ggtitle("(nonAS, IDR<0.1)Ridge Reg(weighted) : median trial PR curve")
  ggsave(paste(outdir,"ridgeReg_ROC_median.pdf",sep=""))
  
  worst.trial <- ridgeStats.df[which(ridgeStats.df$AUC == min(ridgeStats.df$AUC)),"trial"]
  df.tpr.worst <- ridgeTPR.df[which(ridgeTPR.df$trial == worst.trial),]
  ggplot(df.tpr.worst, aes(FPR,TPR))+geom_abline(slope=1) + geom_smooth() + 
    theme_bw()+xlim(0,1)+ylim(0,1)+geom_point()+
    ggtitle("(nonAS, IDR<0.1)Ridge Reg(weighted) : worst trial PR curve")
  ggsave(paste(outdir,"ridgeReg_ROC_worst.pdf",sep=""))
   
}

# ridgeTPR.df[which(ridgeTPR.df$trial == 1),] -> df
getElkansC <- function(probs,label){
  n <- length(probs)
  e1 <- (1/n) * sum(probs[which(label == 1)])
  e2 <- sum(probs[which(label == 1)])/sum(probs)
  e3 <- max(probs)
  data.frame(e1,e2,e3)
}

#Jan 10, 2013 
analyzeRidgeRegressionPU <- function(datadir= getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/"),
                                     outdir  = getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/PULearning/")){
 
  ridgeTPR.df <- readInTable(file=paste(datadir,"ridgeWeightStats-ROC.tab",sep=""))
  ridgeTPR.df$probs <- ridgeTPR.df$dist
  if (!file.exists(outdir)){dir.create(outdir)}
  
  cEstimates.df <- ddply(ridgeTPR.df,.(trial),summarise, 
                      e1 = (1/n) * sum(probs[which(label == 1)]), 
                      e2 = sum(probs[which(label == 1)])/sum(probs), 
                      e3 = max(probs),
                      trial=mean(trial))
  exportAsTable(df=cEstimates.df, file=paste(outdir,"cEstimatesAllTrials.tab",sep=""))
  
  cEstimates.melt <- melt(cEstimates.df,id.var="trial")
  
  ggplot(cEstimates.melt[which(cEstimates.melt$variable != "e3"),], aes(x=value,fill=variable,alpha=0.4)) + geom_density() +
    theme_bw() + xlim(0,max(c(cEstimates.df$e1,cEstimates.df$e2 ))) + 
    geom_vline(xintercept=mean(cEstimates.df$e1),color="red") +
    geom_vline(xintercept=mean(cEstimates.df$e2),color="green") +
    ggtitle("(nonAS, IDR<0.1)Ridge Regression(weighted): 200 trials, \nElkan's c estimate(e1 & e2)")
  ggsave(paste(outdir,"elkansCestimate-Density.pdf",sep=""))
  
  ggplot(ridgeTPR.df,aes(probs))+geom_density() +
    ggtitle("(nonAS, IDR<0.1)Ridge Regression(weighted): 200 trials") +xlim(0,1) +theme_bw()
  ggsave(paste(outdir,"ridgeReg-prob-spectrum.pdf",sep=""))
  
  ddply(ridgeTPR.df, .(trial), transform, 
        probsE1 = probs / ((1/n) * sum(probs[which(label == 1)])),
        probsE1 = probs / (sum(probs[which(label == 1)])/sum(probs))
        )
  
}

#Jan 10, 2013 
ridgeRegCompareEigenVals <- function(datadir= getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/"),
         outdir =  getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/eigenRankCompare/")){
  if (!file.exists(outdir)){dir.create(outdir)}
  
  ## get all the eigen stufff..
  
  df = read.table(getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/fullExpr.tab"),stringsAsFactors=FALSE,header=TRUE)
  doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  
  nolab <- df[which(df$label == 0),]
  lab <- df[which(df$label == 1),]
 eigen.df <-  getEigenVectorsDensity(lab=lab,nolab,cols=exprCols.lpa)
  
  ## get the Ridge Regression data...
  
  ridgeGenePred.df <- readInTable(file=paste(datadir,"genePredictionCounts.tab",sep=""))
  ridgeGenePred.df$lncRnaName <- as.vector(sapply(ridgeGenePred.df[, "gene_id"],
                   function(x)unlist(strsplit(x, "\\."))[1]))
  
  comb.df <- merge(ridgeGenePred.df, eigen.df, by="lncRnaName")
  exportAsTable(df=comb.df, file=paste(outdir,"eigenRankWridgeRegCounts.tab",sep=""))
  
  ggplot(comb.df, aes(rank,count)) + geom_point() + 
    theme_bw() + ggtitle("(nonAS, IDR<0.1)ridge reg vs. eigen val\ncount predicted 1 vs. eigenRank")
  ggsave(paste(outdir,"ridgeEigenCompare-rank-count.pdf",sep=""))
  
  ggplot(comb.df, aes(rank,fill=factor(label))) + geom_density(alpha=I(0.4)) +
    theme_bw() + ggtitle("(nonAS, IDR<0.1)ridge reg vs. eigen val\neigenRank")
  ggsave(paste(outdir,"ridgeEigenCompare-rank-spectrum.pdf",sep=""))
  
  comb.df$weightedCount = ifelse(comb.df$label == 1,
                                 (comb.df$count +1)/sum(comb.df[which(comb.df$label == 1),"count"] +1),
                                 (comb.df$count +1)/sum(comb.df[which(comb.df$label == 0),"count"] +1) )
  
  ggplot(comb.df, aes(rank,weight=weightedCount,fill=factor(label))) + 
    geom_density(alpha=I(0.4)) + 
    theme_bw() + ggtitle("(nonAS, IDR<0.1)ridge reg vs. eigen val\nrank weighted by count(pseudo count of one added)")
  ggsave(paste(outdir,"ridgeEigenCompare-rank-spectrumWeighted.pdf",sep=""))
  
  ridgeTPR.df <- readInTable(file=paste(datadir,"ridgeWeightStats-ROC.tab",sep=""))
  ridgeTPR.df$lncRnaName <- as.vector(sapply(ridgeTPR.df$lncRnaName,
                                                      function(x)unlist(strsplit(x, "\\."))[1]))
  ridgeProbs.df <- ddply(ridgeTPR.df, .(lncRnaName),summarise, aveProbs = mean(probs))
  combProbs.df <- merge(ridgeProbs.df, eigen.df, by="lncRnaName")
  
  ggplot(combProbs.df, aes(x=aveProbs, y = rank,size=label,color=factor(label))) + geom_point() +
    scale_size(range = c(3,5)) + theme_bw() + 
    layer(data=combProbs.df[which(combProbs.df$label == 1),], 
          mapping=aes(x=aveProbs, y = rank, size=4), 
          geom="point") +
    xlab("average probability(ridge reg)") + ylab("eigenRank") +
    ggtitle("(nonAS, IDR<0.1)ridge reg probs vs. eigen val\nrank weighted by count\n(pseudo count of one added)")
  ggsave(paste(outdir,"ridgeEigenCompare-probsVeigenRank.pdf",sep=""))
  
  exportAsTable(df=combProbs.df, file=paste(outdir,"eigenRankWridgeRegProbs.tab",sep=""))
}



doit <- function(){
  rm(list=unlist(ls()))
  source('~/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/fullAnalysis-Jan2014.R')
#   runCompleteAnalysis(outdirFull =  getFullPath("plots/fullAnalysisExperiment/"),
#                        lncGeneIdShorts.vec = c( "IDRlessthan0_1","allLncRNA"),
#                        biotypes.vec = c("antisense","lincRNA","processed_transcript","all_biotypes","lincProcTrans"),
#                        cols.vec = c("lpa","bothPullDowns"))
  runCompleteAnalysis(outdirFull =  getFullPath("plots/fullAnalysisExperiment/"),
                      lncGeneIdShorts.vec = c( "allLncRNA"),
                      biotypes.vec = c("lincProcTrans"),
                      cols.vec = c("bothPullDowns"))
  
  
  
  set.seed(1)
  runLogRegTest()
  
  analyzeLogRegTest()
  
  runRidgeRegression()
  analyzeRidgeReg()
  #getElkansC()
  analyzeRidgeRegressionPU()
  ridgeRegCompareEigenVals()
  testRegByPackage()
}




