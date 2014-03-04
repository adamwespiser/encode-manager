source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R",sep=""))

# http://big.crg.cat/computational_biology_of_rna_processing/lncrna_data

analyzePredOne <- function(outdir= getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/")){

derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                     function(x)unlist(strsplit(x, "\\."))[1]))
lnc.out.file = getFullPath("data/lncExprWithStats_allTrans_idr.tab")
idr.df <- readInTable(lnc.out.file)

outdir 
found.df <- read.csv(paste(outdir,"foundGeneByRidgeRegression.tab",sep=""), stringsAsFactors=FALSE, sep = "\t")

pcFound.df <- found.df[which(found.df$gene_biotype == "protein_coding"),]
pcFound_geneBiotype.df <- pcFound.df[c("ensembl_gene_id", "gene_biotype")]
colnames(pcFound_geneBiotype.df) <- c("ensembl_gene_id", "jan2014_biotype")


derrienPcFound.df <- derrien.df[which(derrien.df$ensembl_gene_id %in% pcFound.df$ensembl_gene_id),]
derrien_geneBiotype.df <- unique(derrienPcFound.df[c("ensembl_gene_id", "gene_biotype")])
colnames(derrien_geneBiotype.df ) <- c("ensembl_gene_id", "derrien_gene_biotype")

merge.df <- merge(derrien_geneBiotype.df, pcFound_geneBiotype.df, by = "ensembl_gene_id")

getExternalGeneIdFromEnsGene(derrienPcFound.df$ensembl_gene_id, 
                             filter="ensembl_gene_id",attr=c("ensembl_gene_id","hgnc_symbol","gene_biotype"))


bm.file= getFullPath("data/lnc_biotype.tab")
bm.df <- read.csv(file=bm.file,sep="\t", stringsAsFactors = FALSE)
bmFound.df <- bm.df[which(bm.df$gene_id_short %in% pcFound.df$ensembl_gene_id),]
bmFound_geneBiotype.df <-  unique(bmFound.df[c("gene_id_short", "bm_biotype", "bm_type")])

bmFound.df$jan_biotype <- getExternalGeneIdFromEnsGene(bmFound.df$gene_id_short, 
                                                       filter="ensembl_gene_id",attr=c("ensembl_gene_id","gene_biotype"))
combined.df <- merge(merge.df,bmFound_geneBiotype.df, by.x= "ensembl_gene_id",by.y="gene_id_short" )
getExternalGeneIdFromEnsGene(combined.df$ensembl_gene_id, 
                             filter="ensembl_gene_id",attr=c("gene_biotype"))

bm_updated.df <- merge(getExternalGeneIdFromEnsGene(bm.df$gene_id_short, 
                             filter="ensembl_gene_id",attr=c("ensembl_gene_id","gene_biotype")),
      bm.df, by.x="ensembl_gene_id",by.y="gene_id_short")


foundWchangedBiotype <- sum(found.df$ensembl_gene_id %in% 
        bm_updated.df[
          which(bm_updated.df$bm_biotype != bm_updated.df$gene_biotype),"ensembl_gene_id"]) / dim(found.df)[1]

foundGeneTypes.df <- data.frame(table(bm_updated.df$gene_biotype))

}

plotRidgeReg <- function(outdir= getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/")){
  
  
  
  
}


createBMfile <- function(){
  
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                                 function(x)unlist(strsplit(x, "\\."))[1]))

  ensembl = useMart("ensembl_mart_46", dataset="hsapiens_gene_ensembl", archive = TRUE)
}


getEnsArch <- function(ens.vec,plotResults=FALSE,
                       outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/")){
  if (missing(ens.vec)){
  
  found.df <- read.csv(paste(outdir,"foundGeneByRidgeRegression.tab",sep=""), stringsAsFactors=FALSE, sep = "\t")
  ens.pc <- found.df[which(found.df$gene_biotype == "protein_coding"), "ensembl_gene_id"]
  ens.vec <- ens.pc
  }
  
  
  #  http://useast.ensembl.org/info/website/archives/index.html
  ens.archive <- list(Ensembl73="sep2013",
                      Ensembl72="jun2013",
                      Ensembl71="apr2013",
                      Ensembl70="jan2013",
                      Ensembl69="oct2012",
                      Ensembl68="jul2012",
                      Ensembl67="may2012",
                      Ensembl66="feb2012",
                      Ensembl65="dec2011",
                      Ensembl64="sep2011",
                      Ensembl63="jun2011",
                      Ensembl62="apr2011",
                      Ensembl61="feb2011",
                      Ensembl59="aug2010",
                      Ensembl54="may2009")
  
  ens.website <- llply(ens.archive, function(x)paste(x,".archive.ensembl.org",sep=""))
  
  
  
  ensembl73=useMart(host='sep2013.archive.ensembl.org', 
                    biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  #listAttributes(ensembl)
  
  
  filter <- "ensembl_gene_id"
  attr <- c("ensembl_gene_id","gene_biotype")
  genes <- derrien.df$ensembl_gene_id
  
  copies.df <- found.df
  for(i in 1:length(names(ens.website))){
    ensBiotypes.df <- getBM(attributes =  c("ensembl_gene_id","gene_biotype"),
                            filters = "ensembl_gene_id", values = ens.vec,
                            mart = useMart(host=ens.website[[i]], 
                                           biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl'),
                            verbose=FALSE)
    colnames(ensBiotypes.df) <- c("ensembl_gene_id", names(ens.website)[i])
    copies.df <- merge(copies.df, ensBiotypes.df, by = "ensembl_gene_id",all.x=TRUE)
    dim(copies.df)
  }
  
  exportAsTable(copies.df, file =paste(outdir, "/lncFoundEnsInfo_ensFull.tab",sep="") )
  copies.melt <- melt(copies.df[,c(1,3,5:19)],id.var="ensembl_gene_id")
  copies.melt$ensemblRelease <- sub(pattern= "Ensembl", replacement = "", copies.melt$variable)
  #copies.melt$ensemblRelease <- as.numeric(substr(copies.melt$variable,8,10))
  
  if (plotResults){
  ggplot(copies.melt, aes(x=factor(ensemblRelease),y=ensembl_gene_id, fill=value,color="black")) + 
    geom_raster() +   theme_bw()
  ggsave(paste(outdir, "/lncFoundEnsBiotypes.pdf",sep=""),height=13,width=8)
  
  ggplot(copies.melt[which(copies.melt$ensembl_gene_id %in% ens.pc),], aes(x=factor(ensemblRelease),y=ensembl_gene_id, fill=value,color="black")) + 
    geom_raster() +   theme_bw() + ggtitle("Ridge Reg id'd w/ protein coding label")
  ggsave(paste(outdir, "/lncFoundEnsBiotypes2.pdf",sep=""),height=13,width=8)
  }
  
}

analyzeEnsemblArchive <- function(outdir){
  if(missing(outdir)){
    outdir = getFullPath("/analysis/ensemblArchive/")
  }
  if(!file.exists(outdir)){
    dir.create(outdir,recursive=TRUE,showWarnings=TRUE)
  }
  df <- readInTable(getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/lncFoundEnsInfo_ensFull.tab"))
  
}

createBmEnsembl73 <- function(){
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                                 function(x)unlist(strsplit(x, "\\."))[1]))
  bm73.file= getFullPath("data/lnc_biotype_ens73.tab")
  
  
  genes <- derrien.df$ensembl_gene_id[1:100]
  ensembl73=useMart(host='sep2013.archive.ensembl.org', 
                  biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  attr <- c("ensembl_gene_id","gene_biotype")
  fullAttrs <- listAttributes(ensembl73)$name[1:20]
  ensBiotypes.df <- getBM(attributes = attr,
                          filters = "ensembl_gene_id", values = derrien.df$ensembl_gene_id,
                          mart = ensembl73, verbose=FALSE)
  
  exportAsTable(ensBiotypes.df, file=bm73.file)
  
}

lncRNAComparisons <- function(outdir){
  expr.df <- readInTable(getFullPath("plots/fullAnalysisExperiment/test/logReg/paramTests/fullExpr.tab"))
  
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                                 function(x)unlist(strsplit(x, "\\."))[1]))
  derrien.df$gene_biotype <- NULL
  
  ensBiotype.df <- readInTable(file=getFullPath("data/lnc_biotype_ens73.tab"))
  
  found.file <- getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/foundGeneByRidgeRegression.tab")
  found.df <- read.csv(file=found.file, stringsAsFactors=FALSE,sep="\t")
  
  expr.df$label <- 0
  expr.df[which(expr.df$lncRnaName %in% found.df$ensembl_gene_id),]$label <- 1
 
  plotMsg = paste("lncRNA group =","lincProcTrans","::: data cols =",columns,"\n(",sum(expr.df$label), "/", 
                  dim(expr.df)[1], ") = (func/total lncRNA) ::: label = pred by ridgeReg", sep=" ")      
  plotLncRNAComparisons(lncDf=expr.df,outdir=paste(outdir,"/comparison/",sep=""),
                        filebase="ridgeRegNonAS-IDR0_1", cols=getLpaColnames(expr.df,getLpaColnames(expr.df,colnames(expr.df))),
                        titleMsg=plotMsg) #need label == 1
  
  
  
}


analyzeEnsVersionsForGeneList <- function(){
    ensBiotype.df <- readInTable(file=getFullPath("data/lncFoundEnsInfo_ensFull.tab"))
  copies.melt <- melt(ensBiotype.df,id.var="ensembl_gene_id")
  copies.melt$ensemblRelease <- sub(pattern= "Ensembl", replacement = "", copies.melt$variable)
  
  
}




analyzegetVersionsForGeneList <- function(indir=getFullPath("data/")){
  outdir <- getFullPath("plots/fullAnalysisExperiment/ensemblVersion/")
  if (!file.exists(outdir)){
    dir.create(outdir,recursive=TRUE,showWarnings=TRUE)
  }
  
  idr.file <- getFullPath("plots/fullAnalysisExperiment/lincProcTrans/lpa-IDRlessthan0_1/comparison/lpa-IDRlessthan0_1.tab")
  idr.df <- readInTable(idr.file)
  gene.idr <- idr.df$external_gene_id
  
  ensv.file =getFullPath("data/lncFoundEnsInfo_ensFull.tab")
  df <- read.table(ensv.file, sep="\t", stringsAsFactors=FALSE,header=TRUE)
  df <- subset(df,ensembl_gene_id %in% gene.idr)
  d <- df
  ens.order <- d[order(factor(d$Ensembl73)),"ensembl_gene_id" ]
  melt.df <- melt(d, id.var="ensembl_gene_id")
  colnames(melt.df) <- c("ensembl_gene_id","ensemblRelease","value")
  melt.df$ensemblRelease <- sub(pattern= "Ensembl", replacement = "", melt.df$ensemblRelease)
  
  ens.pc <- df[which(df$Ensembl73 == "protein_coding"), "ensembl_gene_id"]
  
  melt.df$value <- levels(melt.df$value)[factor(ifelse(melt.df$value %in% c("lincRNA","processed_transcript","protein_coding","antisense","sense_intronic","pseudogene"), melt.df$value,"other"))]
  ggplot(melt.df, aes(x=factor(ensemblRelease),y=ensembl_gene_id, fill=value,color="black")) + 
    geom_tile() +   theme_bw() + scale_y_discrete(limits=ens.order) +
    theme(axis.text.y=element_text(size=5,face="bold")) +
    ggtitle("Ensembl gene biotype:\nIDR < 0.1 & biotype=lincProcTrans")
  ggsave(paste(outdir, "/lncFoundEnsBiotypes.pdf",sep=""),
         height=50,width=8,limitsize=FALSE)
    
  df$plotGroup = sort(rep_len(1:9,dim(df)[1]))
  cols <- c("ensembl_gene_id",colnames(df)[grep(colnames(df), pattern="Ense")])
  
  for(i in 1:9){
    df.local <- subset(df, plotGroup == i)
    df.local <- df.local[,cols]
    ens.order.local <- df.local[order(factor(df.local$Ensembl73)),"ensembl_gene_id" ]
    
    ens.order.local <- df.local[order(factor(df.local$Ensembl73)),"ensembl_gene_id" ]
    melt.local.df <- melt(df.local[cols], id.var="ensembl_gene_id")
    colnames(melt.local.df) <- c("ensembl_gene_id","ensemblRelease","value")
    melt.local.df$ensemblRelease <- sub(pattern= "Ensembl", replacement = "", melt.local.df$ensemblRelease)
    
    ens.pc <- df[which(df$Ensembl73 == "protein_coding"), "ensembl_gene_id"]
    
    melt.local.df$value <- levels(melt.local.df$value)[factor(ifelse(melt.local.df$value %in% c("lincRNA","processed_transcript","protein_coding","antisense","sense_intronic","pseudogene"), 
                                                         melt.local.df$value,"other"))]
    ggplot(melt.local.df, aes(x=factor(ensemblRelease),y=factor(ensembl_gene_id), fill=value,color="black")) + 
      geom_tile() +   theme_bw() + 
      theme(axis.text.y=element_text(size=5,face="bold"))+
      ggtitle("Ensembl gene biotype:\nIDR < 0.1 & biotype=lincProcTrans")
    ggsave(paste(outdir, "/lncFoundEnsBiotypes-part=",i,".pdf",sep=""),
           height=12,width=8,limitsize=FALSE)
  }
}



getMeltedPhastConsForPred <- function(getMelt=TRUE
                                      aggrFN=mean,
                                      found.file = getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/foundGeneByRidgeRegression.tab"),
                                      derrien.file = getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv") ,
                                      idr.file =  getFullPath("plots/fullAnalysisExperiment/lincProcTrans/lpa-IDRlessthan0_1/comparison/lpa-IDRlessthan0_1.tab")
                                      ){
  ### found file
  found.df <- read.csv(file=found.file, stringsAsFactors=FALSE,sep="\t")

  ### get lncRNA used in training/test ridge regression
  idr.df <- readInTable(idr.file)
  gene.idr <- idr.df$external_gene_id
   
  ### conservation scores for each lncRNA transcript
  derrien.df <- read.csv(file=derrien.file,
                         stringsAsFactors=FALSE)
  derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                                 function(x)unlist(strsplit(x, "\\."))[1]))
  derrien.df$gene_biotype <- NULL
  derrienCon.df <- derrien.df[which(derrien.df$ensembl_gene_id %in% gene.idr),c("ensembl_gene_id",
                                colnames(derrien.df)[grep(x=colnames(derrien.df),pattern="phast")])]
  
  ### aggregate the lncRNA conservation by gene
  conMean.df <- ddply(derrienCon.df, .(ensembl_gene_id),numcolwise(aggrFN))  
  
  if (!getMelt){
    conMean.df
  
  } else {
  conMean.df$label <- factor(ifelse(conMean.df$ensembl_gene_id %in% found.df$ensembl_gene_id, 1, 0))
  con.melt <- melt(conMean.df, id.var=c("ensembl_gene_id", "label"), 
                   measure.var=colnames(derrien.df)[grep(x=colnames(derrien.df),pattern="phast")])  
  con.melt$clade_region <- sub(x=sub(x=as.character(con.melt$variable),pattern="phastcons_",replacement=""),
                               pattern="_score",replacement="")
  con.melt$clade <- factor(unlist(lapply(strsplit(x=con.melt$clade_region, split="_"), function(x)x[1])))
  con.melt$region <- factor(unlist(lapply(strsplit(x=con.melt$clade_region, split="_"), function(x)x[2])))
  con.melt
  }
}


lncRNAConservation <- function(
  outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/conservation/")
){
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  con.melt <- getMeltedPhastConsForPred(getMelt=TRUE)
  ggplot(con.melt, aes(x = value, fill = label)) + geom_density(alpha = I(0.4)) +
    facet_grid(clade ~ region) + theme_bw() +
    ggtitle("Ridge Regression predict y=(0,1)\nLincProcTrans, IDR < 0.1")
  ggsave(paste(outdir,"cladeRegionCons_density.pdf", sep = ""))
  
  ggplot(con.melt, aes(y = value, x = label)) + geom_boxplot() +
    facet_grid(clade ~ region) + theme_bw() +
    ggtitle("Ridge Regression predict y=(0,1)\nLincProcTrans, IDR < 0.1")
  ggsave(paste(outdir,"cladeRegionCons_boxplot.pdf",sep = ""))
  
  ggplot(con.melt, aes(x = log(value), fill = label)) + geom_density(alpha = I(0.4)) +
    facet_grid(clade ~ region) + theme_bw() +
    ggtitle("Ridge Regression predict y=(0,1)\nLincProcTrans, IDR < 0.1")
  ggsave(paste(outdir,"cladeRegionCons_LOG-density.pdf",sep = ""))
  
  ggplot(con.melt, aes(y = log(value), x = label))+geom_boxplot() +
    facet_grid(clade ~ region) + theme_bw() +
    ggtitle("Ridge Regression predict y=(0,1)\nLincProcTrans, IDR < 0.1")
  ggsave(paste(outdir,"cladeRegionCons_LOG-boxplot.pdf", sep = ""))
}


