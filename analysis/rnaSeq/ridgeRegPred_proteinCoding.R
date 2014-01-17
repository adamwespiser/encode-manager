source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R",sep=""))

# http://big.crg.cat/computational_biology_of_rna_processing/lncrna_data
derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                     function(x)unlist(strsplit(x, "\\."))[1]))
lnc.out.file = getFullPath("data/lncExprWithStats_allTrans_idr.tab")
idr.df <- readInTable(lnc.out.file)

outdir = getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/")
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


createBMfile <- function(){
  
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  derrien.df$ensembl_gene_id <- as.vector(sapply(derrien.df$LncRNA_GeneId,
                                                 function(x)unlist(strsplit(x, "\\."))[1]))

  ensembl = useMart("ensembl_mart_46", dataset="hsapiens_gene_ensembl", archive = TRUE)
}


testEnsArch <- function(){
  outdir <- getFullPath("plots/fullAnalysisExperiment/test/logReg/ridgeRegression/")
  found.df <- read.csv(paste(outdir,"foundGeneByRidgeRegression.tab",sep=""), stringsAsFactors=FALSE, sep = "\t")
  ens.pc <- found.df[which(found.df$gene_biotype == "protein_coding"), "ensembl_gene_id"]
  
  
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
                            filters = "ensembl_gene_id", values = found.df$ensembl_gene_id,
                            mart = useMart(host=ens.website[[i]], 
                                           biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl'),
                            verbose=FALSE)
    colnames(ensBiotypes.df) <- c("ensembl_gene_id", names(ens.website)[i])
    copies.df <- merge(copies.df, ensBiotypes.df, by = "ensembl_gene_id",all.x=TRUE)
    dim(copies.df)
  }
  
  exportAsTable(copies.df, file =paste(outdir, "/lncFoundEnsInfo.tab",sep="") )
  copies.melt <- melt(copies.df[,c(1,3,5:19)],id.var="ensembl_gene_id")
  copies.melt$ensemblRelease <- sub(pattern= "Ensembl", replacement = "", copies.melt$variable)
  #copies.melt$ensemblRelease <- as.numeric(substr(copies.melt$variable,8,10))
  
  
  ggplot(copies.melt, aes(x=factor(ensemblRelease),y=ensembl_gene_id, fill=value,color="black")) + 
    geom_raster() +   theme_bw()
  ggsave(paste(outdir, "/lncFoundEnsBiotypes.pdf",sep=""),height=13,width=6)
  
  ggplot(copies.melt[which(copies.melt$ensembl_gene_id %in% ens.pc),], aes(x=factor(ensemblRelease),y=ensembl_gene_id, fill=value,color="black")) + 
    geom_raster() +   theme_bw()
  ggsave(paste(outdir, "/lncFoundEnsBiotypes.pdf",sep=""),height=13,width=6)
  
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












