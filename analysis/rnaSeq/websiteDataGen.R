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


getTransDataWithGeneId <- function(){
  # creates the transExpressionForLncRna.tab file, which has transcript expr and the related gene information
expr.data <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/combinedExprWithStats.tab"

derrien.data <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/derrienTransAndGene.tab"
outFile <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/transExpressionForLncRna.tab"

expr.df <- readInTable(expr.data)
read.csv(file=derrien.data,sep="\t",stringsAsFactors=FALSE) -> derrien.df

expr.df[which(expr.df$transcript_id %in% derrien.df$LncRNA_Txid),]-> exprLnc.df
merge(exprLnc.df,derrien.df,by.x="transcript_id",by.y="LncRNA_Txid") -> comb.df

exprCol.names <- read.csv(file=  "/Users/adam/work/research/researchProjects/encode/encode-manager/data/combinedExprWithStats_transEachSample.tab",sep='\t')

comb.out <- comb.df[c("LncRNA_GeneId",colnames(comb.df)[-77])]
exportAsTable(df=comb.out,file=outFile)


}