home <- Sys.getenv("HOME")
projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()


source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R")

#out files created: 


# UNCOMMENT FOR GENE = TOPT transcript for each sample
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
#table(expr.in.df$transcriptType)

#process lncRNAs and export to table
lnc.expr.df <- expr.in.df[which(expr.in.df$transcriptType == "derrienLncRNA"),]
lnc.expr.df <- merge(lnc.expr.df,derrien.df[c("LncRNA_Txid","LncRNA_GeneId")],by.x="transcript_id",by.y="LncRNA_Txid")

#get the rank of transcripts by total expression (sum over all cell-types for transcript)
getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"transcript_id"]),function(x)x[1]), ~ .id )
#ldply(split(lnc.expr.df,lnc.expr.df$LncRNA_GeneId),getTranscriptForMaxCols) -> lnc.expr.maxTransGeneCol.df

# havaPC
cd.expr.df <- expr.in.df[which(expr.in.df$transcriptType == "havanaPC"),]
cd.expr.df <- merge(cd.expr.df,geneMap.df[c("gene_id","transcript_id")])


#getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"transcript_id"]),function(x)x[1]), ~ .id )
#ldply(split(cd.expr.df,cd.expr.df$gene_id),getTranscriptForMaxCols) -> cd.expr.maxTransGeneCol.df

colnames(lnc.expr.df) <- colnames(cd.expr.df)


#read in the condensed cshl rna seq data...
expr.cols <- colnames(expr.in.df[,2:33])
#expr.uniq.cols <- unique(unlist(lapply(colnames(expr.in.df)[2:33],function(x)strsplit(x,".long"))))
#expr.cells = expr.uniq.cols[grep("Poly",expr.uniq.cols,invert=TRUE)]

#get colwise means for lnc and cd...
lnc.gene.max.df <- ddply(lnc.expr.df,.(gene_id),numcolwise(max))
cd.gene.max.df <- ddply(cd.expr.df,.(gene_id),numcolwise(max))


expr.df <- appendAnalysis( rbind(lnc.gene.max.df,cd.gene.max.df), expr.cols)

#test conditions
# expr.tmp <- expr.df
# expr.df <- head(expr.df,20)


#not sure what to do about merging...
#merged.df<-merge(lnc.df,df.2,by.x="LncRNA_Txid",by.y="transcript_id")
  # STEP 2: subset into lncRNA, havana_pc

expr.df$geneType <- "unknown"
expr.df[which(expr.df$gene_id %in% derrien.df$LncRNA_GeneId),]$geneType   <- "derrienLncRNA"
expr.df[which(expr.df$gene_id %in% geneMap.df$gene_id),]$geneType <- "havanaPC"

# -> lengths check out -> all found...

exportAsTable(expr.df,combined.out.file)

#process lncRNAs and export to table
lnc.expr.df <- expr.df[which(expr.df$geneType == "derrienLncRNA"),]
#lnc.expr.df <- merge(lnc.expr.df,derrien.df,by.x="transcript_id",by.y="LncRNA_Txid")

#get the rank of transcripts by total expression (sum over all cell-types for transcript)
#ldply(split(lnc.expr.df[c("transcript_id","sumExpr","LncRNA_GeneId")],lnc.expr.df$LncRNA_GeneId),rankSumExpr)->lnc.expr.rank.df
#lnc.expr.df <- merge(lnc.expr.df,lnc.expr.rank.df[c("transcript_id","rank")])

#get the total number of transcripts per gene:
#ldply(split(lnc.expr.df[c("transcript_id","sumExpr","LncRNA_GeneId")],lnc.expr.df$LncRNA_GeneId),function(df)dim(df)[1]) -> lnc.expr.iso.df
#colnames(lnc.expr.iso.df) <- c("LncRNA_GeneId","isoforms")
#lnc.expr.df <- merge(lnc.expr.df,lnc.expr.iso.df)

# seperate genome coords for lnc
#lnc.coords.df<-ldply(as.list(lnc.expr.df$Hg19_coordinates),function(line)unlist(strsplit(line,'[:-]')))
#colnames(lnc.coords.df)<-c("chr","start","stop")
#lnc.expr.df <- cbind(lnc.expr.df,lnc.coords.df)

exportAsTable(lnc.expr.df,lnc.out.file)

#process coding transcripts and export to table
cd.expr.df <- expr.df[which(expr.df$geneType == "havanaPC"),]
#cd.expr.df <- merge(cd.expr.df,geneMap.df)

#get the rank of each transcript for genes by sum of expression
#ldply(split(cd.expr.df[c("transcript_id","sumExpr","gene_id")],cd.expr.df$gene_id),rankSumExpr)->cd.expr.rank.df
#cd.expr.df <- merge(cd.expr.df,cd.expr.rank.df[c("transcript_id","rank")])

#get the number of isoforms per gene
#ldply(split(cd.expr.df[c("transcript_id","sumExpr","gene_id")],cd.expr.df$gene_id),function(df)dim(df)[1]) -> cd.expr.iso.df
#colnames(cd.expr.iso.df) <- c("gene_id","isoforms")
#cd.expr.df <- merge(cd.expr.df,cd.expr.iso.

exportAsTable(cd.expr.df,cd.out.file)
#iso.df<- data.frame(lnc=c(as.vector(table(lnc.expr.df$rank)),rep(0,24)),pc=as.vector(table(cd.expr.df$rank)),isoforms=seq(1,63))
#melt(iso.df,id.var="isoforms")->iso.melt.df


#exportAsTable(merged.df,getFullPath("/Users/adam/Desktop/Gencode_lncRNAsv7_summaryTable_05_02_2012_lncExression.tab"))

#geneIsoform.df$isoformSqr<-apply(geneIsoform.df[2],1,function(x){x*x})
#geneI.df <- ddply(geneIsoform.df,.(isoforms),subset,sum)
#ggplot(geneIsoform.df,aes(x=isoformSqr))+geom_bar()



###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

test <- function () {
  cd.out.file = getFullPath("data/cdExprWithStats_transEachSample.tab")
  lnc.out.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
  combined.out.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  lnc.expr.df <- data.frame()
  cd.expr.df <- data.frame()
  expr.df <- data.frame()
  expr.file = "/Volumes/MyBook/Volume1/scratch/encodeRnaSeq/cshl/allCellsCombined.space"
  expr.in.df <- read.table(file=expr.file,  header=TRUE, stringsAsFactors=FALSE)
  geneMap.file = "/Volumes/MyBook/Volume1/scratch/gencode/gencode.v7.HAVANA.pc_transcripts.annotation.gtf"
  geneMap.df <- read.table(file=geneMap.file,  header=FALSE, stringsAsFactors=FALSE)
  colnames(geneMap.df) <- c("chr","start", "stop", "transcript_id", "gene_id", "strand")
  derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  expr.in.df$transcriptType <- "unknown"
  expr.in.df[which(expr.in.df$transcript_id %in% derrien.df$LncRNA_Txid),]$transcriptType   <- "derrienLncRNA"
  expr.in.df[which(expr.in.df$transcript_id %in% geneMap.df$transcript_id),]$transcriptType <- "havanaPC"
  lnc.expr.df <- expr.in.df[which(expr.in.df$transcriptType == "derrienLncRNA"),]
  lnc.expr.df <- merge(lnc.expr.df,derrien.df[c("LncRNA_Txid","LncRNA_GeneId")],by.x="transcript_id",by.y="LncRNA_Txid")
  getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"transcript_id"]),function(x)x[1]), ~ .id )
  #ldply(split(lnc.expr.df,lnc.expr.df$LncRNA_GeneId),getTranscriptForMaxCols) -> lnc.expr.maxTransGeneCol.d
  cd.expr.df <- expr.in.df[which(expr.in.df$transcriptType == "havanaPC"),]
  cd.expr.df <- merge(cd.expr.df,geneMap.df[c("gene_id","transcript_id")])
  colnames(lnc.expr.df) <- colnames(cd.expr.df)

  
  lnc.gene.max.df <- ddply(lnc.expr.df,.(gene_id),numcolwise(max))
  cd.gene.max.df <- ddply(cd.expr.df,.(gene_id),numcolwise(max))
  
  expr.df <- appendAnalysis( rbind(lnc.gene.max.df,cd.gene.max.df), colnames(expr.in.df[,2:33]))
 
  exportAsTable(expr.df[which(expr.df$gene_id %in% derrien.df$LncRNA_GeneId),],lnc.out.file)
  exportAsTable(expr.df[which(expr.df$gene_id %in% geneMap.df$gene_id),],cd.out.file)
}
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

createAllTransFile <- function () {
  cd.out.file = getFullPath("data/cdExprWithStats_allTrans.tab")
  lnc.out.file = getFullPath("data/lncExprWithStats_allTrans.tab")
  combined.out.file = getFullPath("data/combinedExprWithStats_allTrans.tab")
  lnc.expr.df <- data.frame()
  cd.expr.df <- data.frame()
  expr.df <- data.frame()
  expr.file = "/Volumes/MyBook/Volume1/scratch/encodeRnaSeq/cshl/allCellsCombined.space"
  expr.in.df <- read.table(file=expr.file,  header=TRUE, stringsAsFactors=FALSE)
  geneMap.file = "/Volumes/MyBook/Volume1/scratch/gencode/gencode.v7.HAVANA.pc_transcripts.annotation.gtf"
  geneMap.df <- read.table(file=geneMap.file,  header=FALSE, stringsAsFactors=FALSE)
  colnames(geneMap.df) <- c("chr","start", "stop", "transcript_id", "gene_id", "strand")
  derrien.file <- getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
  derrien.df <- read.csv(file=getFullPath("data/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv"), stringsAsFactors=FALSE)
  expr.in.df$transcriptType <- "unknown"
  expr.in.df[which(expr.in.df$transcript_id %in% derrien.df$LncRNA_Txid),]$transcriptType   <- "derrienLncRNA"
  expr.in.df[which(expr.in.df$transcript_id %in% geneMap.df$transcript_id),]$transcriptType <- "havanaPC"
  lnc.expr.df <- expr.in.df[which(expr.in.df$transcriptType == "derrienLncRNA"),]
  lnc.expr.df <- merge(lnc.expr.df,derrien.df[c("LncRNA_Txid","LncRNA_GeneId")],by.x="transcript_id",by.y="LncRNA_Txid")
  getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"transcript_id"]),function(x)x[1]), ~ .id )
  #ldply(split(lnc.expr.df,lnc.expr.df$LncRNA_GeneId),getTranscriptForMaxCols) -> lnc.expr.maxTransGeneCol.d
  cd.expr.df <- expr.in.df[which(expr.in.df$transcriptType == "havanaPC"),]
  cd.expr.df <- merge(cd.expr.df,geneMap.df[c("gene_id","transcript_id")])
  colnames(lnc.expr.df) <- colnames(cd.expr.df)
  
  
  #lnc.gene.max.df <- ddply(lnc.expr.df,.(gene_id),numcolwise(max))
  #cd.gene.max.df <- ddply(cd.expr.df,.(gene_id),numcolwise(max))
  
  expr.df <- appendAnalysis( rbind(lnc.expr.df,cd.expr.df), colnames(expr.in.df[,2:33]))
  
  exportAsTable(expr.df[which(expr.df$gene_id %in% derrien.df$LncRNA_GeneId),],lnc.out.file)
  exportAsTable(expr.df[which(expr.df$gene_id %in% geneMap.df$gene_id),],cd.out.file)
}
