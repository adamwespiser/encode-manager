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
source( paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/exprLib.R",sep=""))
source( paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R",sep=""))


sourceFile <- getFullPath("data/lncExprWithStats_allTrans_idr.tab")
df <- readInTable(sourceFile)
df <- df[which(df$averageExpr != 0),]
df[["gene_id_short"]] <- sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})


createDir = function(outdir){
if(!file.exists(outdir)){
  dir.create(outdir)
}}

removeExprColTags = function(x){
  sub(pattern="(\\.longPolyA|\\.longNonPolyA)\\.(RPKM2|RPKM1|IDR|COMB)",replacement="",x)
}


makeTranscriptPerGene = function(){
  #remove "."
  sourceFile <- getFullPath("data/lncExprWithStats_allTrans_idr.tab")
  df <- readInTable(sourceFile)
  expr.cols <- colnames(df)[grep(x=colnames(df),pattern="Poly")]
  cell.types <- unique(sapply(expr.cols,removeExprColTags))
  
  
  df.melt <- melt(df, by=c("transcript","gene_id","transcriptType"))
  df.melt$celltype <- sapply(df.melt$variable,function(x)unlist(strsplit(as.character(x),"\\."))[1])
  df.melt$pulldown <- apply(df.melt$variable,function(x)unlist(strsplit(as.character(x),"\\."))[2])
  df.melt$measure <- sapply(df.melt$variable,function(x)unlist(strsplit(as.character(x),"\\."))[3])
  df.melt$transCelltype <- apply(df.melt,1,function(x)paste(x[["transcript"]],x[["celltype"]],x[["pulldown"]],sep=""))
  df.melt$geneCelltype <- apply(df.melt,1,function(x)paste(x[["gene_id"]],x[["celltype"]],x[["pulldown"]],sep=""))
  df.melt.comb <- df.melt[which(df.melt$measure == "COMB"),]
  df.melt.comb.transPerGene <- ddply(df.melt.comb,.(geneCelltype),function(df)df[which(df$value == max(df$value))[1],])
  #df.melt.comb.transPerGene = ddply(df.melt.comb,.(celltype,gene_id,),subset,value == max(value))
  #df.melt.comb.transPerGene = ldply(split(df.melt.comb,df.melt.comb$geneCelltype),function(df)df[which(df$value == max(df$value)),])
  df.melt.topTrans <- df.melt[which(df.melt$transCelltype %in% df.melt.comb.transPerGene$transCelltype),]
  exportAsTable(df.melt.topTrans,file="./data/lncRnaExpr_2reps_idr_melt.tab")
 
  
  df.melt.topTrans$one <- 1
  df.topTrans <- dcast(df.melt.topTrans[c("one","gene_id","variable","value")], gene_id ~ variable,sum)  
  
  exportAsTable(df.topTrans,file="./data/lncRnaExpr_2reps_idr.tab")
  
  df.topTrans.split <- dcast(df.melt.topTrans[c("measure","celltype","pulldown","gene_id","value")], gene_id + celltype + pulldown ~  measure,sum) 
  exportAsTable(df.topTrans.split,file="./data/lncRnaExpr_2reps_idr_split.tab")
  #test 
  d$geneCelltype <- apply(d,1,function(x)paste(x[["gene_id"]],x[["celltype"]],sep=""))
  
  table(df$gene_id)[which(table(df$gene_id) > 5) ]
  ensg <- "ENSG00000255248.1"
  df.topTrans[which(df.topTrans$gene == ensg),] -> t
  df[which(df$gene == ensg),] -> tt
  
  
  
  
  exps <- as.vector(sapply(colnames(df)[grep("COMB",colnames(df))], function(x)gsub(x=x,pattern="\\.COMB",replacement="")))
  getCols <- function(x)paste0(x,c(".COMB",".RPKM1",".RPKM2","IDR"))
  cellTypeList <- llply(exps,function(ex)df.topTrans[,c("gene",getCols(ex))])
  
  
}


preProcessData = function(file.split=getFullPath("data/lncRnaExpr_2reps_idr_split.tab"),
                       file.long = getFullPath("/data/lncRnaExpr_2reps_idr.tab")){
  
  df <- readInTable(file.split)
  celltypes <- unique(df$celltype)  
  df[["gene_id_short"]] <-  sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  

  df.long <- readInTable(file.long)
  df.long[["averageExpr"]] <- apply(df.long[,colnames(df.long)[grep(x=colnames(df.long),pattern="COMB")]],1,function(x)mean(x))
  df.long[["gene_id_short"]] <-  sapply(as.character(df.long$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  
 # df.long <- df.long[which(df.long$averageExpr != 0),]
  ensg.expressed <- df.long$gene_id_short  
  
  ens.df <- getEnslist()
  ensg.functional <- unique(ens.df$ensembl_gene_id)
  ensg.extern.df <- unique(ens.df[c("ensembl_gene_id","external_gene_id")])
  
  df <- df[which(df$gene_id_short %in% ensg.expressed),]
  df$label <- rep(0,length(df$gene_id_short))
  df[df$gene_id_short %in% ensg.functional,"label"] <- 1 
  
  # add external id 
  df.one <- df[which(df$label == 1),]
  df.zero <- df[which(df$label == 0),]
  df.zero$external_gene_id <- df.zero$gene_id_short
  df.one <- merge(df.one,ensg.extern.df, by.x = "gene_id_short", by.y = "ensembl_gene_id")
  df.comb <- rbind(df.one[,colnames(df.one)],df.zero[,colnames(df.one)])
  df.comb <- within(df.comb,{
    
    exprNone = (RPKM1 == 0 & RPKM2 == 0)
    exprOne  = (RPKM1 > 0 & RPKM1 == 0) | (RPKM2 > 0 & RPKM1 == 0)
    exprBoth = (RPKM1 > 0 & RPKM2 > 0)
    RPKMSum = RPKM1 + RPKM2
  })
  df.comb
  
}

# ds
getLabelFn <- function(cutoff){
  function(df, celltype, pulldown ){
    df.tmp = df[which(df$exprOne == TRUE & df$celltype == celltype & df$pulldown == pulldown),]
    cutoff = quantile(df.tmp[["RPKMSum"]],c(cutoff))[[1]]
    sum(df[which(df$celltype == celltype & df$pulldown == pulldown & df$COMB > cutoff ),"label"])
  }
}


labelFnIDR <- function(cutoff,df, celltype, pulldown ){
    sum(df[which(df$celltype == celltype & df$pulldown == pulldown & df$IDR <= cutoff ),"label"])
}



# normalize two RPKMs
normRPKM_first <- function(r1, r2){
  (((r1 + r2)/2) * 10^3)/ (max(r1,r2)/2)
}


normRPKM <- function(RPKM1, RPKM2){
(((RPKM1 + RPKM2)/2) *1000)/max(RPKM1,RPKM2)/2
}

plotCummFnByVariable <- function(outdir = getFullPath("plots/rnaSeq-QC/allCells/"),
                                 df = preProcessData()){
  
  
  
}


plotData <- function(outdir = getFullPath("plots/rnaSeq-QC/allCells/")){
  
  if(!file.exists(outdir)){
    dir.create(outdir)
  }
  
  df = preProcessData()
  celltypes = unique(df$celltype)
  percentiles <- c(0.5,0.7,0.8,0.9,0.95)
  df.out <- expand.grid(c("longPolyA","longNonPolyA"),celltypes)
  colnames(df.out) <- c("pulldown", "celltype")
  for (p in c(0,0.1,0.25,0.5,0.8,0.9,0.95,0.99)){
    f = getLabelFn(p)
    df.out[[gsub(x=paste("p",p,sep=""),replacement="_",pattern="\\.")]] <- apply(df.out,1,function(x)f(df=df,celltype=x[["celltype"]],pulldown=x[["pulldown"]]))
  }
  
  exportAsTable(df.out,file = paste(outdir,"lncRNA_foundInEitherCutoff_appliedToAllLncRNA.tab",sep=""))
  
  ########### IDR
  df.out <- expand.grid(c("longPolyA","longNonPolyA"),celltypes)
  colnames(df.out) <- c("pulldown", "celltype")
  
  
  for (p in c(0.0,0.01,0.04,0.08,0.1,0.15,0.20,0.25,0.3,0.4,0.5,0.7,1,10)){
    df.out[[gsub(x=paste("p",p,sep=""),replacement="_",pattern="\\.")]] <- apply(df.out,1,function(x)labelFnIDR(cutoff = p,df=df,celltype=x[["celltype"]],pulldown=x[["pulldown"]]))
  }
  
  exportAsTable(df.out,file = paste(outdir,"lncRNA_IDRcutoff_appliedToAllLncRNA.tab"),sep="")
  
  
  
  
  outdir.rpkm <- paste(outdir,"rpkmCompare/",sep="")
  createDir(outdir.rpkm)
  for (cell in celltypes){
    outfile <- paste(outdir.rpkm,cell,"_PolyA_RPKM1_2.pdf")
    df.local <- df[which(df$celltype == cell & df$pulldown == "longPolyA"),]
    plot.title = paste("RPMKs for celltype=",cell, "\npulldown=longPolyA")
    
    g <- ggplot(df.local[order(df.local$label),], aes(x = RPKM1, y = RPKM2, color = factor(label), size= label )) + 
      geom_point() +
      theme_bw() + scale_size(range=c(2,5))
    g + ggtitle(plot.title)
    ggsave(paste(outdir.rpkm,cell,"_PolyA_RPKM1_2.pdf"),height=6,width=8)
    # zoom in plots 
    g + xlim(0,quantile(df$RPKM1,0.95)) + 
      ylim(0,quantile(df$RPKM2,0.95)) +
      ggtitle(paste(plot.title,"\ntop 5% of datapoints removed"))
    ggsave(paste(outdir.rpkm,cell,"_PolyA_RPKM1_2_zoom.pdf",sep=""),height=6,width=8)
    
    ggplot(df.local, aes(x=IDR, fill = factor(label))) + geom_bar() + 
      theme_bw() + 
      geom_vline(x=0.1) + 
      facet_wrap(~label,ncol=1,scale="free_y")+
      ggtitle(paste(cell,"longPolyA\nIDR (label: 1=function,0=unlab)"))
    ggsave(paste(outdir.rpkm,cell,"_PolyA_IDR.pdf",sep=""),height=6,width=8)
    
  }
  
  ddply(df, .(celltype), summarize, exprBoth=sum(exprBoth),exprOne=sum(exprOne), exprNone=sum(exprNone))
  
  
  #v as.data.frame(table(df[df$IDR < 0.4, "celltype"]))[["Freq"]]
  
  ggplot(df,aes(x=IDR,fill=factor(label)))+geom_density(alpha=I(0.4)) + theme_bw() +
    ggtitle("All celltype/lnc IDR entries")
  ggsave(paste(outdir,"IDR_density.pdf",sep=""))
  
  ggplot(df, aes(x=IDR,fill=factor(label)))+geom_bar(binwidth=0.02)+theme_bw()+
    facet_wrap(~label,ncol=1,scale="free_y")+
    ggtitle("IDR distribution of lncRNA s.t. IDR > 0\nfacet 0 = unlabel & facet 1 = label")
  ggsave(paste(outdir,"IDR_distro.pdf",sep=""))
  
  
  
  
  
  df.lpa <- df[which(df$pulldown == "longPolyA"),]
  idr.df = as.data.frame(table(df.lpa[df.lpa$IDR < 1, "celltype"]))
  idr.df$idr = 1
  for(i in 2:1000){
    p = (i - 1) / 1000
    tmp = as.data.frame(table(df.lpa[df.lpa$IDR < p, "celltype"]))
    tmp$idr = p
    idr.df = rbind(idr.df, tmp)
  }
  
  dflabel = df.lpa[which(df.lpa$label == 1),]
  idrlabel.df = as.data.frame(table(dflabel[dflabel$IDR < 1, "celltype"]))
  idrlabel.df$idr = 1
  for(i in 2:1000){
    p = (i - 1)/ 1000
    tmp = as.data.frame(table(dflabel[dflabel$IDR < p,"celltype"]))
    tmp$idr = p
    idrlabel.df = rbind(idrlabel.df, tmp)
  }
  
  
  idrlabel.df$label = "label"
  idr.df$label = "unlabel"
  idrratio.df = idrlabel.df
  idrratio.df$Freq = idrlabel.df$Freq / idr.df$Freq
  idrratio.df$label = "labelOverUnlabel"
  idr.comb = rbind(idr.df,idrlabel.df,idrratio.df)
  
  colnames(idr.comb) <- c("celltype","lncRnaFound","idrCutoff","label")
  
  ggplot(idr.comb, aes(x=idrCutoff,y=lncRnaFound,color=celltype))+geom_line()+
    theme_bw()+xlim(0,0.5)+
    facet_wrap(~label,ncol=1,scale="free_y")
  ggsave(paste(outdir,"idrCuttof_plusRatio.pdf",sep=""))
  
   #################### COMB
  comb.df = as.data.frame(table(df.lpa[df.lpa$COMB > 0, "celltype"]))
  comb.df$comb = 0
  for(i in 2:300){
    p = (i - 1) / 100
    tmp = as.data.frame(table(df.lpa[df.lpa$COMB > p, "celltype"]))
    tmp$comb = p
    comb.df = rbind(comb.df, tmp)
  }
  
  dflabel = df[which(df.lpa$label == 1),]
  comblabel.df = as.data.frame(table(dflabel[dflabel$COMB > 0, "celltype"]))
  comblabel.df$comb = 0
  for(i in 2:300){
    p = (i - 1)/ 100
    tmp = as.data.frame(table(dflabel[dflabel$COMB > p,"celltype"]))
    tmp$comb = p
    comblabel.df = rbind(comblabel.df, tmp)
  }
  
  
  comblabel.df$label = "label"
  comb.df$label = "unlabel"
  combratio.df = comblabel.df
  combratio.df$Freq = comblabel.df$Freq / comb.df$Freq
  combratio.df$label = "labelOverUnlabel"
  together.df = rbind(comb.df,comblabel.df,combratio.df)
  
  colnames(together.df) <- c("celltype","lncRnaFound","combinedRPKMCutoff","label")
  
  ggplot(together.df, aes(x=combinedRPKMCutoff,y=lncRnaFound,color=celltype))+geom_line()+
    theme_bw()+xlim(0,3)+
    facet_wrap(~label,ncol=1,scale="free_y") +
    ggtitle("Combined RPKM distribution of lncRNA s.t. IDR > 0\nfacet 0 = unlabel & facet 1 = label")
  ggsave(paste(outdir,"combinedRPKM_cutoffPlusRatio.pdf",sep=""))
  
 
  # number of lncRNA found in each cell type
  c.df <- as.data.frame(table(df[which(df$COMB > 0), "celltype"]))
  c.df$src = "combinedRPKM"
  l.df <- as.data.frame(table(df[which(df$IDR <  0.5), "celltype"]))
  l.df$src = "IDR"
  cl.df <- rbind(c.df,l.df)
  colnames(cl.df) <- c("celltype","foundLncRNA", "measure")
  ggplot(cl.df, aes(x=celltype, y= foundLncRNA,fill=measure)) + 
    geom_histogram(position="dodge",stat="identity") +
    ggtitle("Number of lncRNA found with COMB > 0, IDR < 0.5")
  ggsave(paste(outdir,"lncRNAFoundByCutoff.pdf",sep=""))
  
  
  
  c.df <- as.data.frame(table(df[which(df$COMB > 0), "celltype"]))
  c.df$src = "combinedRPKM"
  l.df <- as.data.frame(table(df[which(df$IDR <  0.1), "celltype"]))
  l.df$src = "IDR"
  cl.df <- rbind(c.df,l.df)
  colnames(cl.df) <- c("celltype","foundLncRNA", "measure")
  ggplot(cl.df, aes(x=celltype, y= foundLncRNA,fill=measure)) + 
    geom_histogram(position="dodge",stat="identity") +
    ggtitle("Number of lncRNA in celltype w/\n RPKM COMB > 0 & IDR < 0.1") + theme_bw() + coord_flip()
  ggsave(paste(outdir,"lncRNAFoundByCutoff-0_1.pdf",sep=""))
}




