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


sourceFile = getFullPath("data/lncExprWithStats_allTrans_idr.tab")
df = readInTable(sourceFile)
df = df[which(df$averageExpr != 0),]
df[["gene_id_short"]] =  sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})





removeExprColTags = function(x){
  sub(pattern="(\\.longPolyA|\\.longNonPolyA)\\.(RPKM2|RPKM1|IDR|COMB)",replacement="",x)
}


makeTranscriptPerGene = function(){
  #remove "."
  sourceFile = getFullPath("data/lncExprWithStats_allTrans_idr.tab")
  df = readInTable(sourceFile)
  expr.cols = colnames(df)[grep(x=colnames(df),pattern="Poly")]
  cell.types = unique(sapply(expr.cols,removeExprColTags))
  
  
  df.melt = melt(df, by=c("transcript","gene_id","transcriptType"))
  df.melt$celltype =sapply(df.melt$variable,function(x)unlist(strsplit(as.character(x),"\\."))[1])
  df.melt$pulldown =sapply(df.melt$variable,function(x)unlist(strsplit(as.character(x),"\\."))[2])
  df.melt$measure = sapply(df.melt$variable,function(x)unlist(strsplit(as.character(x),"\\."))[3])
  df.melt$transCelltype = apply(df.melt,1,function(x)paste(x[["transcript"]],x[["celltype"]],sep=""))
  df.melt.comb  = df.melt[which(df.melt$measure == "COMB"),]
  df.melt.comb.transPerGene = ddply(df.melt.comb,.(celltype,gene_id),subset,value == max(value))
  df.melt.topTrans = df.melt[which(df.melt$transCelltype %in% df.melt.comb.transPerGene$transCelltype),]
  exportAsTable(df.melt.topTrans,file="./data/lncRnaExpr_2reps_idr_melt.tab")
  df.topTrans = dcast(df.melt.topTrans[c("gene_id","transcript","variable","value")], transcript + gene_id ~ variable)
  exportAsTable(df.topTrans,file="./data/lncRnaExpr_2reps_idr.tab")
  
    

  
  
  
}