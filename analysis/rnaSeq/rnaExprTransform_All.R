home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)

source(getFullPath("analysis/RnaSeq/exprLib.R"))

#out files created: 
cd.out.file = getFullPath("data/cdExprWithStats.tab")
lnc.out.file = getFullPath("data/lncExprWithStats.tab")
combined.out.file = getFullPath("data/combinedExprWithStats.tab")

#gencode.v7.HAVANA.pc_transcripts.annotation.gtf
#gencode.file <- "/Users/adam/work/research/researchProjects/encode/encode-manager/data/g"


#file with combined expression 
expr.file = "/Volumes/MyBook/Volume1/scratch/encodeRnaSeq/cshl/allCellsCombined.space"
expr.df <- read.table(file=expr.file,  header=TRUE, stringsAsFactors=FALSE)

# get a list of data.frame(transcript=... , gene= ...) for determining canonical transcript note: somewhere in data file on cluster?
geneMap.file = "/Volumes/MyBook/Volume1/scratch/gencode/gencode.v7.HAVANA.pc_transcripts.annotation.gtf"
geneMap.df <- read.table(file=geneMap.file,  header=FALSE, stringsAsFactors=FALSE)
colnames(geneMap.df) <- c("chr","start", "stop", "transcript_id", "gene_id", "strand")

# get derrien's list of lncRNA
derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)



#read in the condensed cshl rna seq data...
expr.cols <- colnames(expr.df[,-1])
expr.uniq.cols <- unique(unlist(lapply(colnames(expr.df)[-1],function(x)strsplit(x,".long"))))
expr.cells = expr.uniq.cols[grep("Poly",expr.uniq.cols,invert=TRUE)]

#test conditions
# expr.tmp <- expr.df
# expr.df <- head(expr.df,20)
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
expr.df <- assignExprLabels(expr.df,cols=expr.cols)

expr.df$expressed = "no"
expr.df[which(expr.df$maxExprType != "none"),]$expressed <- "yes"

#not sure what to do about merging...
#merged.df<-merge(lnc.df,df.2,by.x="LncRNA_Txid",by.y="transcript_id")
  # STEP 2: subset into lncRNA, havana_pc

expr.df$transcriptType <- "unknown"
expr.df[which(expr.df$transcript_id %in% derrien.df$LncRNA_Txid),]$transcriptType   <- "derrienLncRNA"
expr.df[which(expr.df$transcript_id %in% geneMap.df$transcript_id),]$transcriptType <- "havanaPC"

# -> lengths check out -> all found...

exportAsTable(expr.df,combined.out.file)

#process lncRNAs and export to table
lnc.expr.df <- expr.df[which(expr.df$transcriptType == "derrienLncRNA"),]
lnc.expr.df <- merge(lnc.expr.df,derrien.df,by.x="transcript_id",by.y="LncRNA_Txid")

#get the rank of transcripts by total expression (sum over all cell-types for transcript)
ldply(split(lnc.expr.df[c("transcript_id","sumExpr","LncRNA_GeneId")],lnc.expr.df$LncRNA_GeneId),rankSumExpr)->lnc.expr.rank.df
lnc.expr.df <- merge(lnc.expr.df,lnc.expr.rank.df[c("transcript_id","rank")])

#get the total number of transcripts per gene:
ldply(split(lnc.expr.df[c("transcript_id","sumExpr","LncRNA_GeneId")],lnc.expr.df$LncRNA_GeneId),function(df)dim(df)[1]) -> lnc.expr.iso.df
colnames(lnc.expr.iso.df) <- c("LncRNA_GeneId","isoforms")
lnc.expr.df <- merge(lnc.expr.df,lnc.expr.iso.df)

# seperate genome coords for lnc
lnc.coords.df<-ldply(as.list(lnc.expr.df$Hg19_coordinates),function(line)unlist(strsplit(line,'[:-]')))
colnames(lnc.coords.df)<-c("chr","start","stop")
lnc.expr.df <- cbind(lnc.expr.df,lnc.coords.df)
exportAsTable(lnc.expr.df,lnc.out.file)

#process coding transcripts and export to table
cd.expr.df <- expr.df[which(expr.df$transcriptType == "havanaPC"),]
cd.expr.df <- merge(cd.expr.df,geneMap.df)

#get the rank of each transcript for genes by sum of expression
ldply(split(cd.expr.df[c("transcript_id","sumExpr","gene_id")],cd.expr.df$gene_id),rankSumExpr)->cd.expr.rank.df
cd.expr.df <- merge(cd.expr.df,cd.expr.rank.df[c("transcript_id","rank")])

#get the number of isoforms per gene
ldply(split(cd.expr.df[c("transcript_id","sumExpr","gene_id")],cd.expr.df$gene_id),function(df)dim(df)[1]) -> cd.expr.iso.df
colnames(cd.expr.iso.df) <- c("gene_id","isoforms")
cd.expr.df <- merge(cd.expr.df,cd.expr.iso.

exportAsTable(cd.expr.df,cd.out.file)



iso.df<- data.frame(lnc=c(as.vector(table(lnc.expr.df$rank)),rep(0,24)),pc=as.vector(table(cd.expr.df$rank)),isoforms=seq(1,63))
melt(iso.df,id.var="isoforms")->iso.melt.df


#exportAsTable(merged.df,getFullPath("/Users/adam/Desktop/Gencode_lncRNAsv7_summaryTable_05_02_2012_lncExression.tab"))

#geneIsoform.df$isoformSqr<-apply(geneIsoform.df[2],1,function(x){x*x})
#geneI.df <- ddply(geneIsoform.df,.(isoforms),subset,sum)
#ggplot(geneIsoform.df,aes(x=isoformSqr))+geom_bar()
