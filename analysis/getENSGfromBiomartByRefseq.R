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
library(biomaRt) # http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
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
getExternalGeneIdFromEnsGene <- function(gene_id,filter="ensembl_gene_id",attr=c("ensembl_gene_id","external_gene_id")){
  getBM(attributes = attr,
               filters = filter, values = gene_id,
               mart = mart)
}

getLncRNADiseasdbVec <- function(file=getFullPath("data/lncRnaAndDisease_refSeq.txt")){
  df <- read.table(file=file,stringsAsFactors=FALSE,header=FALSE)
  colnames(df) <- c("refseq")
  sapply(df$refseq,function(x)getEnsGFromRef(value=splitOnDot(x),mart=mart))->f.vec
  unique(as.data.frame(sapply(f.vec,function(x)x[1]))[,1])
}

getLncRNAdbInteract <- function(file=getFullPath("data/lncRNAdbInteract_REFSEQ.tab")){

 df <- read.csv(file=file,sep="\t",stringsAsFactors=FALSE,header=FALSE)  
 colnames(df) <- c("refseq","lncRnaName")
 f.vec = sapply(df$refseq,function(x)getEnsGFromRef(value=splitOnDot(x),mart=mart))
 f.out = as.character(unique(as.data.frame(sapply(f.vec,function(x)x[1]))[,1]))
 f.out = f.out[!is.na(f.out)]
 f.out
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

createSourcesForLnc <- function(in.f=getFullPath("data/lncList_ensemblIdsFoundLncRna2.list"),
                                out.f=getFullPath("data/lncEvidenseOfFunction.tab")){
  disease <-     getLncRNADiseasdbVec()
  lncRnaDb <-    getLncRnaDbFromExport()
  lncList <-     getlncdListDf()
  genenames <-   getEnsFromGeneNameList()
  funLncDb <-    getLncsGenesFromFunctionalLncDb()
  lncDbScrape <- getLncsFromLncRnaDbWebScrap()
  lit2 <-        getLncsFromLitSearch2()
  interact <-    getLncRNAdbInteract()
  all.found = c(lit2[["ensembl_gene_id"]],lncList[["gene_id"]],interact,lncDbScrape[["ensembl_gene_id"]],lncRnaDb,funLncDb[["ensembl_gene_id"]],as.character(disease),genenames)
 
  
  lnc.df <- readInTable(file=in.f)
  
  
  evidense.df <- data.frame( ensembl_gene_id = unique(all.found))
  derrien.df = read.csv(getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"))
  #annot.df = getAnnotLncDf()
  evidense.df = within(evidense.df,{
   
    litSearch = (ensembl_gene_id %in% lit2[["ensembl_gene_id"]] | ensembl_gene_id %in% lncList[["gene_id"]])
    lncInteract =  ensembl_gene_id %in% interact
    lncRnaDb = (ensembl_gene_id %in% lncDbScrape[["ensembl_gene_id"]] | ensembl_gene_id %in% lncRnaDb)
    lncFuncDb = ensembl_gene_id %in% funLncDb[["ensembl_gene_id"]]
    lncDiseaseDb = ensembl_gene_id %in% as.character(disease)
    lncGeneNamesDb = ensembl_gene_id %in% genenames
    foundInDerrien = ensembl_gene_id %in% sapply(as.character(derrien.df$LncRNA_GeneId),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
    
    sum = ifelse(foundInDerrien,1,0)*100+ ifelse(lncRnaDb,1,0)*1.4 + ifelse(lncFuncDb,1,0)*1.3 + ifelse(lncDiseaseDb,1,0)*1.2 + ifelse(lncGeneNamesDb,1,0) * 1.1 + ifelse(litSearch,1,0)*1.1 + ifelse(lncInteract,1,0)*1.05
    
  })
  
  unique(lnc.df$ensembl_gene_id)
  
  evidense.df[order(evidense.df$sum),"rank"] = seq_along(evidense.df$sum)

  g <- getEnslist()
  g = g[which(g$ensembl_gene_id %in% evidense.df$ensembl_gene_id),]
  ename.df <- merge(evidense.df,unique(data.frame(name=g$external_gene_id,ensembl_gene_id=g$ensembl_gene_id)),by="ensembl_gene_id")
  exportAsTable(ename.df,file=out.f)
  
  #ename.df <- transform( ename.df,
  #                       ensembl_gene_id = ordered(ensembl_gene_id, levels = sort(-table(sum))))
  
  #ename.df <- arrange(ename.df,sum)
  ename.df$ensembl_gene_id <- factor(ename.df$ensembl_gene_id , levels=arrange(ename.df,sum)[["ensembl_gene_id"]], ordered=TRUE)
  ename.df$name <- factor(ename.df$name , levels=arrange(ename.df,sum)[["name"]], ordered=TRUE)
  ename.melt = melt(ename.df,id.var=c("rank","ensembl_gene_id","name","sum"))
  #ename.melt$ensembl_gene_id <- factor(ename.df$ensembl_gene_id , levels=arrange(ename.df,sum)[["ensembl_gene_id"]], ordered=TRUE)
  ggplot(melt(ename.df,id.var=c("rank","ensembl_gene_id","name","sum")),aes(y=name,x=reorder(variable,sum),fill=ifelse(value,0,1))) + 
    geom_tile()+scale_fill_gradient(low = "#00002A",high = "#E0F5FF") + 
    theme(axis.text.y=element_text(size=3,color="black"),axis.text.x=element_text(size=3,color="black"))
  
  ggsave(file="/Users/adam/work/research/researchProjects/encode/encode-manager/plots/lncCompare/lncRnaFunctional/sourceOfLncs.pdf",,height=7,width=5)
  
  lncDb.count1 =  dim(read.csv(file=getFullPath("data/lncRNAdbNrCodesHuman.list"),stringsAsFactors=FALSE,header=FALSE))[1]
  
  lncLit.count = lncDb.count1 + dim(read.table(file="/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncdList_parsed.tab",stringsAsFactors=FALSE,header=TRUE))[1]
  lncDisease.count =dim(read.table(file= "/Users/adam/work/research/researchProjects/encode/encode-manager/data/lncRnaAndDisease_LNC-REFSEQ.tab",stringsAsFactors=FALSE,header=FALSE))[1]
  lncDis.count =  dim(read.table(file=getFullPath("data/lncRnaAndDisease_refSeq.txt"),stringsAsFactors=FALSE,header=FALSE))[1]
  
  
  
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
  

getEnsemblBiotypeDerrienList <- function(){
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                                 function(x)unlist(strsplit(x, "\\."))[1]))
  #listMarts(archive=TRUE)
  ensembl73=useMart(host='sep2013.archive.ensembl.org', 
                    biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  #listAttributes(ensembl)
  
  
  filter <- "ensembl_gene_id"
  attr <- c("ensembl_gene_id","gene_biotype")
  genes <- derrien.df$ensembl_gene_id
  ensBiotypes.df <- getBM(attributes =  c("ensembl_gene_id","gene_biotype"),
                          filters = "ensembl_gene_id", values = genes,
                          mart = ensembl73,
                          verbose=TRUE)
  ensBiotypes.df
  
}



















  

