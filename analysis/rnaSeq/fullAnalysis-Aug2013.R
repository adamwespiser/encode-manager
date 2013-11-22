#The purpose of this script is to run: lncRNA comparison, PCA, logistic Regression, and eigenRank 
# the output will be directed to one directory
#
#
#
home <- Sys.getenv("HOME")

source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/eigenRank_08272013.R",sep=""))
source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/rnaSeq_qc.R",sep=""))
source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/logisticRegPcaExprData.R",sep=""))
#source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R",sep=""))

#LOG Reg 

editStatsForLncDf <- function(expr.df, cols){
  
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,entropy,"entropyExpr")
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,sum,"sumExpr")
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,var,"varianceExpr")
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,mean,"averageExpr")
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,min,"minExpr")
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,max,"maxExpr")
  
  expr.df <- applyAndAppendFnDf(expr.df,expr.cols,function(x)maxExpr(x,expr.cols),"maxExprType")
  expr.df <- applyAndAppendFnDf(df=expr.df,cols=expr.cols,FUN=function(x)threeTimesRestSpecifity(x,expr.cols),"threeTimesRest")
  expr.df <- applyAndAppendFnDf(df=expr.df,cols=expr.cols,FUN=function(x)nTimesRestSpecifity(x,expr.cols),"nTimesRest")
  expr.df <- applyAndAppendFnDf(df=expr.df,cols=expr.cols,FUN=function(x)calcTissSpec(x,expr.cols),"tissSpec")
  expr.df
  
}



plotLncRNAComparisons <- function(lncDf=local.df,outdir=outdir,
                                   cols=expr.cols ,titleMsg="",
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
  melt.df <- melt(lncDf[c("gene_id",expr.cols,id.vec)], 
                  measure.vars=sort(expr.cols),
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

runCompleteAnalysis <-function(outdir =  getFullPath("plots/fullAnalysisExperiment/")){
  
  #get biotype info
  bm.file= getFullPath("data/lnc_biotype.tab")
  bm.df <- read.csv(file=bm.file,sep="\t", stringsAsFactors = FALSE)
  biotypes.list = list(antisense = bm.df[which(bm.df$bm_biotype == "antisense"),"gene_id_short"],
                       processed_transcript = bm.df[which(bm.df$bm_biotype == "processed_transcript"),"gene_id_short"],
                       lincRNA = bm.df[which(bm.df$bm_biotype == "lincRNA"),"gene_id_short"],
                       all_biotypes = bm.df[["gene_id_short"]],
                       remove_antisense =  bm.df[which(bm.df$bm_biotype != "antisense"),"gene_id_short"]
  )
  
  biotypes.vec = c("antisense","lincRNA","processed_transcript","all_biotypes","remove_antisense")
  
  
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
  
  #for(columns in c("lpa","lnpa","bothPullDowns")){
  
  for (biotype in biotypes.vec){
    lncRNA.biotype <- biotypes.list[[biotype]]
    basedir = paste(outdir, "/", biotype,"/", sep = "")
    for(columns in c("lpa","bothPullDowns","lnpa")){ 
      expr.cols <- cols.list[[columns]]
      for(lncGroup in lncGeneIdShorts.vec){
        
        
        
        ## EigenRank
        expr.rows.group <- rows.list[[lncGroup]]
        expr.rows <- expr.rows.group[which(expr.rows.group %in% lncRNA.biotype )]
        
        expr.cols <- cols.list[[columns]]
        print(paste("starting",columns,lncGroup))
        outdir = paste(basedir,columns,"-",lncGroup,sep="")
        
        if(!file.exists(outdir)){dir.create(outdir)}
        outdir = paste(outdir,"/",sep="")
        
        
        local.df = df[which(df$gene_id_short %in% as.vector(expr.rows)),]
        local.df = local.df[which(!apply(local.df[expr.cols],1,function(x)sum(x)) == 0),]
        
        
        plotMsg = paste("lncRNA group =",lncGroup,"::: data cols =",columns,"\n(",sum(local.df$label), "/", 
                        dim(local.df)[1], ") = (func/total lncRNA) ::: biotype =",biotype, sep=" ")      
        #    plotLncRNAComparisons <- function(lncDf=local.df,outdir=outdir,
        #                                      filename=paste(columns,biotype,sep="-"),cols=expr.cols ,titleMsg="",
        #                                      foundColword="",filebase=""){
        
        if (identical(analysisPref$compare,TRUE)){
          print("comparison") 
          plotLncRNAComparisons(lncDf=local.df,outdir=paste(outdir,"/comparison/",sep=""),
                                filebase=paste(columns,lncGroup,sep="-"), cols=expr.cols,
                                titleMsg=plotMsg)
        }
        
        
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
        ## Log. Regression
        if (identical(analysisPref$logreg,TRUE)){
          print("logistic regression")
          
          bestCols.df <- data.frame(
            cols = exprCols,
            score = sapply(expr.cols,function(x){mean(local.df[which(local.df$label == 1),x]) /mean(local.df[which(local.df$label == 0),x])}))
          topNCols = 10
          exprColsBestIndex = bestCols.df[order(bestCols.df$score,decreasing=TRUE)[1:topNCols],"cols"]
          exprLogReg = expr.cols[which(exprCols %in% exprColsBestIndex)] 
          
          
          # logisticRegPcaExprData.R
          runLogReg(lncDf=local.df,outdir=paste(outdir,"logReg",sep=""),
                    cols=exprLogReg,iter=5,debug= FALSE,
                    filebase=paste(columns,lncGroup,sep="-"),
                    titleMsg=plotMsg)
        } # end log reg
      }
    }
    
    
  }
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















