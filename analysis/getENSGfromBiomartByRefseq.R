home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

#source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R")
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)

library(plyr)
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
#                  filters = "refseq_ncrna", values = "NR_001463.4",
#                  mart = mart)


splitOnDot <- function(string){
  
  as.vector(strsplit(string,"\\.")[[1]])
}

getEnsFromRef <- function(value,filter="refseq_ncrna",mart=mart){
  # must have mart object init'd ###### mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl ")###### 
  getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
        filters = filter, values = value,
        mart = mart)
  }


getEnsGFromRef <- function(value,filter="refseq_ncrna",mart=mart){
  # must have mart object init'd ###### mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl ")###### 
  
  getBM(attributes = c("ensembl_gene_id"),
        filters = filter, values = value,
        mart = mart)[1]
}

getLncRNADiseasdbVec <- function(file=getFullPath("data/lncRnaAndDisease_refSeq.txt")){
  df <- read.table(file=file,stringsAsFactors=FALSE,header=FALSE)
  colnames(df) <- c("refseq")
  sapply(df$refseq,function(x)getEnsGFromRef(value=splitOnDot(x),mart=mart))->f.vec
  unique(as.data.frame(sapply(f.vec,function(x)x[1]))[,1])
}

getLncRNAdbDf <- function(file=getFullPath("data/lncRnaAndDisease_refSeq.txt")){
  
  
} 
  
getlncDiseaseDf <- function(f = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncRnaAndDisease_LNC-REFSEQ.tab"){
  df <- read.table(file=f,stringsAsFactors=FALSE,header=FALSE)
  colnames(df)  <- c("lncRnaName","refseq")
  df$gene_id <- apply(df,1,function(x)getEnsGFromRef(x[2],mart=mart))
  df[!is.na(df$gene_id),]
}

getlncdListDf <- function(f = "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncdList_parsed.tab"){
  df <- read.table(file=f,stringsAsFactors=FALSE,header=TRUE)
  colnames(df) <- c("lncRnaName","gene_id","foundInDerrien")
  df
}


getOrthosDf <- function(f="/Users/adam/work/research/researchProjects/encode/encode-manager/data/mouse.all.rep20.ENST.txt"){
  df <- read.table(file=f,stringsAsFactors=FALSE,header=TRUE)
  colnames(df) <- "transcript_id"
  derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  data.frame(lncRnaName="mouse?",gene_id=derrien.df[derrien.df$LncRNA_Txid %in% df$transcript_id,"LncRNA_GeneId"],foundInDerrien="true")
  
}  
  
getDerriensAntiSenseDf <- function(f=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")){
  derrien.df <- read.csv(file=f, stringsAsFactors=FALSE)
  df <- derrien.df[derrien.df[["gene_biotype"]] == "antisense",]
  data.frame(
            lncRnaName =  df$Hg19_coordinates,
            gene_id = df$LncRNA_GeneId,
            foundInDerrien= TRUE
            )
  
}

getLncRnaDbFromExport <- function(f=getFullPath("data/lncRNAdbNrCodesHuman.list")){
  df <- read.csv(file=f,stringsAsFactors=FALSE,header=FALSE)
  r.df<-ldply(df$V1,function(nd_code)getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","refseq_ncrna","description","external_gene_id"),filters= c("refseq_ncrna"),values=list(nd_code),mart=mart))
  unique(r.df$ensembl_gene_id)
}


getEnsFromGeneNameList <- function(genenameFile=getFullPath("data/geneNamesLNCRNA.list")){
  df <- getEnsFromGeneNameListFull()
  df[["ensembl_gene_id"]]
}

getLncsGenesFromFunctionalLncDb <- function(file=getFullPath("./data/funcLncRnaDb_scrape.list")){
  funcLncE.list <- ldply(readInTable(getFullPath("/data/funcLncRnaDb_scrape.list"))[["entry"]],
                         function(nd_code)getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","description","external_gene_id","hgnc_symbol"),
                                                filters= c("hgnc_symbol"),
                                                values=list(nd_code),
                                                mart=mart))
  funcLnc.list <- ldply(readInTable(getFullPath("/data/funcLncRnaDb_scrape.list"))[["entry"]],
                        function(nd_code)getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","description","external_gene_id","hgnc_symbol"),
                                               filters= c("wikigene_name"),
                                               values=list(nd_code),
                                               mart=mart))
  
  combined.df <- unique(rbind(funcLnc.list,funcLncE.list))
  combined.df <- combined.df[which(combined.df$gene_biotype != "protein_coding"),]
  combined.df
  
}

getLncsFromLncRnaDbWebScrap <- function(file=getFullPath("/data/lncRnaDbWebScrapForCoords.list")){
  df <- read.csv(file,stringsAsFactor=FALSE)
  lncDbScrape.df <- ldply(df[["coord"]],
                        function(nd_code)getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","description","external_gene_id","hgnc_symbol"),
                                               filters= c("chromosomal_region"),
                                               values=list(unlist(strsplit(nd_code,"chr"))[2]),
                                               mart=mart))
  
  
  
  lncDbScrape.df[which(lncDbScrape.df$gene_biotype != "snoRNA" & lncDbScrape.df$gene_biotype != "protein_coding"),] -> lncDbScrape.df
  
}

getLncsFromLitSearch2 <- function(file=getFullPath("data/lncRnaLitSearch.tab")){
 # df <- read.csv(file,stringsAsFactor=FALSE)
  lit2 <- ldply(readInTable(getFullPath("data/lncRnaLitSearch.tab"))[["gene_id"]],
                         
                function(nd_code)getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","description","external_gene_id","hgnc_symbol"),
                                                filters= c("ensembl_gene_id"),
                                                values=list(nd_code),
                                                mart=mart))
  
  
}







getEnsFromGeneNameListFull <- function(genenameFile=getFullPath("data/geneNamesLNCRNA.list")){
  df <- readInTable(genenameFile)
  gn.df<-ldply(df[[1]],
              function(nd_code)getBM(attributes = c("ensembl_gene_id","gene_biotype","description","external_gene_id","hgnc_symbol"),
                                     filters= c("refseq_ncrna"),
                                     values=list(nd_code),
                                     mart=mart))
  gn.df
}
createEnsList <- function(f=getFullPath("data/lncList_ensemblIdsFoundLncRna.list")){
  out.f=getFullPath("data/lncList_ensemblIdsFoundLncRna.list")
 
  disease <-     getLncRNADiseasdbVec()
  lncRnaDb <-    getLncRnaDbFromExport()
  lncList <-     getlncdListDf()
  genenames <-   getEnsFromGeneNameList()
  funLncDb <-    getLncsGenesFromFunctionalLncDb()
  lncDbScrape <- getLncsFromLncRnaDbWebScrap()
  lit2 <-        getLncsFromLitSearch2()
  
  all.vec <- unique(c(lit2[["ensembl_gene_id"]], lncDbScrape[["ensembl_gene_id"]],as.character(disease),lncRnaDb,genenames,lncList[["gene_id"]],funLncDb[["ensembl_gene_id"]]))
  lncs <- unique(as.vector(sapply(all.vec,function(s)unlist(strsplit(s,"\\."))[1])))
  r.df<-ldply(lncs[!is.na(lncs)],
              function(nd_code)getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","gene_biotype","description","external_gene_id","hgnc_symbol","refseq_ncrna"),
                                     filters= c("ensembl_gene_id"),
                                     values=list(nd_code),
                                     mart=mart))
  #r.df[which(r.df$gene_biotype != "protein_coding"),]
  exportAsTable(r.df,file=out.f)

}

getEnslist <- function(file=getFullPath("data/lncList_ensemblIdsFoundLncRna2.list")){
  if (!file.exists(file)){
    createEnsList()
  }
  else{
    readInTable(file)
  }
  
}

getEnsListAttrRich <- function(file=getFullPath("data/lncList_ensemblIdsFoundLncRna2.list")){
  ens.df <- getEnslist()
  r.df<-ldply(unique(ens.df$ensembl_gene_id),
              function(nd_code)getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","gene_biotype","description","external_gene_id","hgnc_symbol","refseq_ncrna","go_id","embl","wikigene_id"),
                                     filters= c("ensembl_gene_id"),
                                     values=list(nd_code),
                                     mart=mart))
  unique(r.df)
}



pcaWithSubset <- function(lnc,annot,exprCols=2:33){
  
  lnc.pca         <- princomp(lnc[exprCols])
  lnc.pca.factors <-  as.data.frame(lnc.pca$scores)
  lnc.pca.factors[["gene_id"]] <- lnc[["gene_id"]]
  lnc.pca.factors[["withinSubset"]] <-lnc[["withinSubset"]]
  lnc.pca.factors[["lncRnaName"]] <-lnc[["lncRnaName"]]
  
  lnc.pca.factors[["Comp.1"]] <-lnc.pca.factors[["Comp.1"]] - min(lnc.pca.factors[["Comp.1"]]) + 1
  lnc.pca.factors[["Comp.2"]] <-lnc.pca.factors[["Comp.2"]] - min(lnc.pca.factors[["Comp.2"]]) + 1
  #lnc.pca.melt <- melt(lnc.pca.factors,id.vars=c("gene_id","withinSubset"))
#  ggplot(lnc.pca.factors,aes(x=Comp.1,fill=Comp.2))+geom_density()
 # ggsave("~/Desktop/pcaCompsForLncTopGene.pdf")
  ggplot(lnc.pca.factors,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)
  
  lnc.sub <- which(lnc.pca.factors$withinSubset == "true")
  lnc.pca.factors <- lnc.pca.factors(union(which(lnc.pca.factors[["Comp.1"]] >0), which(lnc.pca.factors[["Comp.2"]]>0)))
  lnc.min <- lnc.pca.factors[which(lnc.pca.factors[["Comp.1"]] != min(lnc.pca.factors[["Comp.1"]])),]

  ggplot(lnc.min,aes(x=log(Comp.1),y=log(Comp.2),color=withinSubset,shape=withinSubset,size=3))+geom_point()+theme_bw() + 
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2)),geom="point",size=4)+
    layer(data=lnc.pca.factors[lnc.pca.factors$withinSubset == "true" ,],mapping=aes(x=log(Comp.1),y=log(Comp.2),hjust=0.9,vjust=-1,label=lncRnaName),geom="text",size=4)
  
  
}






evalEnsList <- function(){
  r.df <- createEnsList()
  cd.out.file = getFullPath("data/cdExprWithStats_transEachSample.tab")
  lnc.out.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
  
  
  # create a list of all transcripts 
  #cd.out.file = getFullPath("data/cdExprWithStats_allTrans.tab")
  #lnc.out.file = getFullPath("data/lncExprWithStats_allTrans.tab")
  
  
  #combined.out.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  #create output for all 
  combined.out.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  
  
  
  #gencode.v7.HAVANA.pc_transcripts.annotation.gtf
  #gencode.file <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/g"
  
  lnc.expr.df <- data.frame()
  cd.expr.df <- data.frame()
  expr.df <- data.frame()
  #file with combined expression 
  expr.file = "/Volumes/MyBook/Volume1/scratch/encodeRnaSeq/cshl/allCellsCombined.space"
  expr.in.df <- read.table(file=expr.file,  header=TRUE, stringsAsFactors=FALSE)
  
  # get a list of data.frame(transcript=... , gene= ...) for determining canonical transcript note: somewhere in data file on cluster?
  geneMap.file = "/Volumes/MyBook/Volume1/scratch/gencode/gencode.v7.HAVANA.pc_transcripts.annotation.gtf"
  geneMap.df <- read.table(file=geneMap.file,  header=FALSE, stringsAsFactors=FALSE)
  colnames(geneMap.df) <- c("chr","start", "stop", "transcript_id", "gene_id", "strand")
  
  # get derrien's list of lncRNA
  derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  
  # figure out transcript membership...
  expr.in.df$transcriptType <- "unknown"
  expr.in.df[which(expr.in.df$transcript_id %in% derrien.df$LncRNA_Txid),]$transcriptType   <- "derrienLncRNA"
  expr.in.df[which(expr.in.df$transcript_id %in% geneMap.df$transcript_id),]$transcriptType <- "havanaPC"
  within(expr.in.df,{trans_short = sapply(transcript_id,function(x)unlist(strsplit(x,"\\."))[[1]])})  -> expr.in.df  #unlist(strsplit(nd_code,"chr"))[2]
}
  
  
  
  

