
#################################
#################################
home <- Sys.getenv("HOME")

projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

#boiler plate helpers

source(getFullPath("analysis/getENSGfromBiomartByRefseq.R"))

garber.base.dir  <- getFullPath("data/garber2014data/")
garber.expr.file <- getFullPath("data/garber2014data/garber2014supp.dat")
garber.analysis.file <- getFullPath("data/garber2014data/garber2014supp-analysis.dat")
garber.ortholog.ens73.file <- getFullPath("data/garber2014data/orthologs-ens73.tab")


getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)
#################################
#################################

source(getFullPath("analysis/rnaSeq/exprLib.R"))


gatherDataFromWeb <- function(datadir = "/home/wespisea/sandbox/gtest/"){
  if (missing(datadir)){
    datadir <- garber.base.dir
  }
  if (!file.exists(datadir)){
    dir.create(datadir,recursive=TRUE)
  }

  print(datadir)
  url.base = "http://garberlab.umassmed.edu/data/humanlincRNAEvolSupplement/"
  url.files.vec <-  c("hominid.expressed.lincs.bed.gz", "mammalian.expressed.bed.gz", "online-supplement.tar.gz")
  url.vec <- paste0(url.base, url.files.vec)
  temp.vec <- sapply(seq_along(url.vec),function(x) tempfile())
  untar
  invisible(sapply(seq_along(url.vec),function(x)download.file(url=url.vec[x],destfile=temp.vec[x])))
  untar(tarfile=temp.vec[3],
        exdir=paste(datadir,"/",sep=""))
  
  invisible(sapply(seq_along(url.vec[1:2]), 
  function(i)write(readLines(con=gzfile(temp.vec[i]),n=-1),
        file=paste(datadir,"/",sub(x=url.files.vec[i],pattern=".gz",replacement=""),sep=""))))
  
  file.copy(from=paste(datadir,"/online-data/supplementary-data.csv",sep=""),
            to=paste(datadir,"/garber2014supp.dat",sep=""))
}

#gatherDataFromWeb(datadir = "/home/wespisea/sandbox/gtest/")


stripGeneVersion <- function(id){
  as.vector(sapply(id,
                 function(x)unlist(strsplit(x, "\\."))[1]))
}

createGarberExpr <- function(){
  df <- read.csv(garber.expr.file, stringsAsFactors=FALSE, sep=";")
  df$ensembl_gene_id <- stripGeneVersion(df$id)
  df.analysis <- editStatsForLncDf(df,getTissueIndex(df))
  exportAsTable(df = df.analysis, file=garber.analysis.file )
  df.analysis
}
#editStatsForLncDf


readInGarberExpr <- function(){
   if (file.exists(garber.analysis.file)){
    df <- read.csv(garber.analysis.file, stringsAsFactors=FALSE, sep="\t")
  }
   else {
     df <- createGarberExpr()
   }
  df
}

getGarberInfo <- function(){
  scan(file=getFullPath("/data/garber2014data/online-data/info.txt"),  character(0),sep="\n",quote=NULL)
}


getTissue <- function(){
  c("brain", "cerebellum", "colon", "heart", "kidney", "liver", "lung", "skm", "spleen", "testes") 
}

getTissueIndex <- function(df){
  cols <- colnames(df)
  which(cols %in% getTissue())
}

getGenomeFromSpecies <- function(species){
  if(missing(species)){
    warning("species is missing -> returning NULL")
  }
  species.list <- list(rat= "rn4", rhesus = "rheMac2",
                       chimp = "panTro3", mouse="mm9",
                       cow="bosTau6")
  
  genome <- species.list[[species]]
  if(is.null(genome)){
    warning("garberData.R::getGenomeFromSpecies -> species not found")
  }
  genome
}

# field for using getBM
getEnsemblOrthologFieldFromaName <- function(species){
  if(missing(species)){
    warning("species is missing -> returning NULL")
  }
  species.list <- list(rat= "rnorvegicus_homolog_ensembl_gene", rhesus = "mmulatta_homolog_ensembl_gene",
                       chimp = "ptroglodytes_homolog_ensembl_gene", mouse="mmusculus_homolog_ensembl_gene",
                       cow="btaurus_homolog_ensembl_gene")
  
  genome <- species.list[[species]]
  if(is.null(genome)){
    warning("garberData.R::getGenomeFromSpecies -> species not found")
  }
  genome
}

getEnsemblOrthologFields <- function(){
c("rnorvegicus_homolog_ensembl_gene", "mmulatta_homolog_ensembl_gene",
  "ptroglodytes_homolog_ensembl_gene", "mmusculus_homolog_ensembl_gene",
  "btaurus_homolog_ensembl_gene")
}

fetchHumanOrthologs <- function(outfile = garber.ortholog.ens73.file){
  g.df <- readInGarberExpr()
  g.human.df <- subset(g.df, species == "hg19")
  ortho.df <- fetchEnsembl73(values=g.human.df$ensembl_gene_id,
                             filter= "ensembl_gene_id",
                             attr=c("ensembl_gene_id",getEnsemblOrthologFields()))
  write.csv(file = garber.ortholog.ens73.file, x = ortho.df)
  ortho.df
}

readInGarberOrtho <- function(){
  if (file.exists(garber.ortholog.ens73.file)){
    df <- read.csv(garber.ortholog.ens73.file, stringsAsFactors=FALSE, sep=",")
  }
  else {
    df <- fetchHumanOrthologs()
  }
  df
}




