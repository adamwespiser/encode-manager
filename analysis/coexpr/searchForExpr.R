library(SRAdb)
#sqlfile <- getSRAdbFile()
#source("http://bioconductor.org/biocLite.R")
#biocLite("SRAdb")
setwd("~/work//research/researchProjects//encode/encode-manager/")
sqlfile <- "/home/wespisea/data/SRAmetadb.sqlite" 
sra_con <- dbConnect(SQLite(),sqlfile)
sra_tables <- dbListTables(sra_con)

hpc.system <- function(cmd){
  r.cmd <- pipe(paste("ssh",ghpc,"\"",cmd,"\""))
  str <- readLines(r.cmd)
  flush(r.cmd)
  close(r.cmd)
  rm(r.cmd)
  str
}

hpc.file.exists <- function(file){
  cmd <- paste("file ",file)
  result <- hpc.system(cmd)
  g.result <- grep(x=result, pattern="No such file or directory")
  if(any(g.result)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}


test <- function(){
  dbListFields(sra_con,"run")
  sqliteQuickSQL(sra_con,"PRAGMA TABLE_INFO(sra)")
  
  
  dbListFields(sra_con,"sra")
  rs <- dbGetQuery(sra_con, paste( "SELECT study_type AS StudyType,
 count( * ) AS Number FROM sra GROUP BY study_type order
 by Number DESC ", sep=""))
  
  rs.species <- dbGetQuery(sra_con, paste( "SELECT taxon_id,common_name,
 count( * ) AS Number FROM sra GROUP BY study_type order
 by Number DESC ", sep=""))
  
  rs <- dbGetQuery(sra_con, paste( "SELECT * FROM sra ", sep=""))
  rs.use <- rs[which(rs$library_strategy == "RNA-Seq"  ),]
  
  
  file <- "./data/sraSearchTerm=Human.csv"
  df <- read.csv(file=file)
  
  df.selection <- df[which(df$Library.Strategy == "RNA-Seq" &
                             df$Library.Source == "TRANSCRIPTOMIC" & 
                             df$Library.Selection == "cDNA" &
                             df$Organism.Name == "Homo sapiens"),]
  
  dim(df.selection)[1]
  
  
  dbGetQuery(sra_con, paste( "SELECT * FROM study", sep=""))      -> study.df
  dbGetQuery(sra_con, paste( "SELECT * FROM experiment", sep="")) -> expr.df
  dbGetQuery(sra_con, paste( "SELECT * FROM sample", sep=""))     -> sample.df
  dbGetQuery(sra_con, paste( "SELECT * FROM run", sep=""))        -> run.df
  dbGetQuery(sra_con, paste( "SELECT * FROM metaInfo", sep=""))   -> metaInfo.df
  dbGetQuery(sra_con, paste( "SELECT * FROM sra", sep=""))        -> sra.df
  dbGetQuery(sra_con, paste( "SELECT * FROM run", sep=""))        -> run.df
  
  
  study.rna.df <- sra.df[which(sra.df$library_strategy == "RNA-Seq" & 
                                 sra.df$taxon_id == 9606),]
  
  
  dbGetQuery(sra_con, paste( "SELECT * FROM experiment where experiment_accession == \"DRX000003\" ", sep=""))
  
  listSRAfile('SRR578654',sra_con)
}


getRNASeqRuns <- function(){
  sqlfile <- "/home/wespisea/data/SRAmetadb.sqlite" 
  sra_con <- dbConnect(SQLite(),sqlfile)
  
  dbGetQuery(sra_con, paste( "SELECT * FROM sra", sep=""))        -> sra.df
  select.df <- sra.df[which(sra.df$library_strategy == "RNA-Seq" & 
                              sra.df$taxon_id == 9606),]
  outfile <- "~/sandbox/sradbExprs.tab"
  write.csv(file=outfile,select.df[exprInfoCols],sep="\t")
  out.short <- "~/sandbox/sradbExprs1.tab"
  shortCols <- c('center_project_name', 'experiment_title','study_accession')
  write.csv(file=out.short,select.df[shortCols],sep="\t")
  
  study.vec <- c("SRP000996","ERP000546", "SRP000941","SRP001371","SRP007461","SRP012682","SRP017465","SRP018838","SRP021478","SRP029899","SRP035524")
  sraC.df <- select.df[which(select.df$study_accession %in% study.vec),]
}

convertToList <- function(sa){
  nms.raw <- as.character(unlist(sapply(unlist(strsplit(sa,split=" \\|\\| ")), function(x)unlist(strsplit(x,split=": "))[1])))
  nms <- as.character(unlist(sapply(nms.raw,function(x)gsub(x=x,pattern=" ",replacement="_"))))
  vals <- as.character(unlist(sapply(unlist(strsplit(sa,split=" \\|\\| ")), function(x)unlist(strsplit(x,split=": "))[2])))
  l <- as.list(vals)
  names(l) <- nms
  l
}

# for sra.df$read_spec
convertToDf <- function(str,recSep,labelSep){
  records <- as.character(unlist(strsplit(str, recSep)))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
  vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
  l <- as.list(vals)
  names(l) <- nms
  as.data.frame(l)
}


sampleAttrToDf <- function(sa.vec){
  Reduce("rbind",sapply(sa.vec,function(x) as.data.frame(convertToList(x))))
}

readInGTExAllMeta <- function(){
  df <- read.csv(file = "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/GTExSraDB-metainfo.tab",stringsAsFactors=FALSE,sep="\t")
  datadir <- "/data/wespisea/gtex/fastq/"
  df$read1 <- paste0(datadir,df$run_accession,"_1.fasta.gz")
  df$read2 <- paste0(datadir,df$run_accession,"_2.fasta.gz")
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$have1 <- file.exists(df$fastq1)
  df$have2 <- file.exists(df$fastq2)
  
  df$cddownloadCmd <- paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df$downloadCmd <- paste0("/home/wespisea/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df$gtexId <- as.character(unlist(sapply(df$sample_attribute, function(x)convertToList(x)[["submitted_sample_id"]])))
  
  df
}

gtex.annot <- "/home/wespisea/data/GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt"
getGTExAnnot <- function(file=gtex.annot){
  df <- read.csv(file,sep="\t",stringsAsFactors=FALSE)
  df$SAMPID <- as.character(sapply(df$SAMPID, function(x)gsub(x=x,pattern="[\\_\\-]",replacement="\\.")))
  df
}

getGTExWithFullAnnot <- function(){
  d <- getGTExAnnot()
  m <- readInGTExAllMeta()
  d$SAMPIDdash = gsub(d$SAMPID, pattern="\\.",replacement="-")
  m$gtexId= sapply(m$gtexId,toupper)
  d$SAMPIDdash= sapply(d$SAMPIDdash,toupper)
  comb <- merge(d,m, by.x="SAMPIDdash",by.y="gtexId")
}



createDownloadFile <- function(){
  df <- readInGTExAllMeta()
  write(paste0(df$downloadCmd,rep(c(" &"," "),length=length(df$downloadCmd))), file="~/sandbox/downloadGTEx.sh") 
}

downloadFileMissing <- function(){
  df <- readInGTExAllMeta()
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$haveFiles <- file.exists(df$fastq1)
  df.need <- df[which(df$haveFiles == FALSE),]
  
  write(paste0(df.need$cddownloadCmd,rep(c(" &"," "),length=length(df.need$cddownloadCmd))), file="~/sandbox/downloadGTEx_need.sh") 
  # http://gap-upload.ncbi.nlm.nih.gov/E2004057-252A-4BC9-ADA3-2463E4AC9B98/SRR612551.sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47
  o <- paste0("wget -O /data/wespisea/gtex/sra/",df.need$run_accession,".sra 'http://gap-upload.ncbi.nlm.nih.gov/E2004057-252A-4BC9-ADA3-2463E4AC9B98/",
              df.need$run_accession,".sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47' --continue",
              rep(c(" &"," "),length=length(df.need$run_accession)))
  
  len <- length(o)
  midpoint <- floor(len/2)
  write(o[1:midpoint], file = "~/sandbox/downloadGTEx_wget1.sh")
  write(o[(midpoint+1):len], file = "~/sandbox/downloadGTEx_wget2.sh")
}

downloadFileMissing_url2 <- function(){
  df <- readInGTExAllMeta()
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$haveFiles <- file.exists(df$fastq1)
  df$haveFiles2 <- file.exists(df$fastq2)
  df.need <- df[which(df$haveFiles == FALSE),]
  
  write(paste0(df.need$downloadCmd,rep(c(" &"," "),length=length(df.need$downloadCmd))), file="~/sandbox/downloadGTEx_need.sh") 
  # http://gap-upload.ncbi.nlm.nih.gov/E2004057-252A-4BC9-ADA3-2463E4AC9B98/SRR612551.sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47
  o <- paste0("wget -O /data/wespisea/gtex/sra/",df.need$run_accession,".sra 'http://gap-upload.ncbi.nlm.nih.gov/D0F68BBC-70E4-4163-BA50-8CFFF55E4CDE/",
              df.need$run_accession,".sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47' --continue",
              rep(c(" &"," "),length=length(df.need$run_accession)))
  
  len <- length(o)
  midpoint <- floor(len/2)
  write(o[1:midpoint], file = "~/sandbox/downloadGTEx2_wget1.sh")
  write(o[(midpoint+1):len], file = "~/sandbox/downloadGTEx2_wget2.sh")
}

getMeanDownloadSpeed <- function(downloadFile,sampleSize=10){
  mean(diff(sapply(1:sampleSize,function(x){Sys.sleep(1);file.info(downloadFile)$size}))/1000)
}
# getMeanDownloadSpeed("/data/wespisea/gtex//fastq/SRR818749_1.fastq.gz")





getGTExSra <- function(){
  sqlfile <- "/home/wespisea/data/SRAmetadb.sqlite" 
  stopifnot(file.exists(sqlfile))
  
  sra_con <- dbConnect(SQLite(),sqlfile)
  
  dbGetQuery(sra_con, paste( "SELECT * FROM sra", sep=""))        -> sra.df
  select.df <- sra.df[which(sra.df$library_strategy == "RNA-Seq" & 
                              sra.df$taxon_id == 9606),]
  
  # get the GTEx samples
  gtex.df <- select.df[which(select.df$study_accession %in% "SRP012682"),]
  # get all the GTEx sample attributes
  sampleAttr   <- sapply(gtex.df$sample_attribute, function(x)convertToList(x))
  # take only SRA GTEx entries with 16 memebers in sample attributes
  gtex.16.df <- gtex.df[which(as.numeric(unlist(sapply(sampleAttr,length))) == 16),]
  gtex.11.df <-  gtex.df[which(as.numeric(unlist(sapply(sampleAttr,length))) == 11),]
  # get sample attributes where sample attrs has 16 members
  # sampleAttr.16   <- sapply(gtex.16.df$sample_attribute[1:100], function(x)as.data.frame(convertToList(x)))
  
  # convert sample attributes to df
  sampleAttr.df <- ldply(gtex.16.df$sample_attribute, function(x)as.data.frame(convertToList(x)))
  
  file.df <- ldply(gtex.16.df$run_accession, function(x)listSRAfile(x,sra_con))
  
  comb.df <- cbind(gtex.16.df,sampleAttr.df,file.df)
}

plotGTEx <- function(){
  comb.df <- getGTExSra()
  outdir <- getFullPath("plots/coexpr/GTEx")
  
  ggplot(comb.df, aes(as.numeric(gap_subject_id),body_site,fill=sex))+geom_tile() +theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb.df$body_site, comb.df$body_site, length))))+
    ggtitle("GTEx data from SRAdb")
  ggsave(paste(outdir,"bodysite-heatmap.pdf",sep=""),height=6,width=12)
  
  ggplot(comb.df, aes(as.numeric(gap_subject_id),histological_type,fill=sex))+geom_tile()+theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb.df$histological_type, comb.df$histological_type, length))))+
    ggtitle("GTEx data from SRAdb")
  ggsave(paste(outdir,"histological_type-heatmap.pdf",sep=""),height=5,width=12)
  
  ggplot(comb.df, aes(histological_type,fill=sex))+geom_bar() + coord_flip() + theme_bw() +
    scale_x_discrete(limits=names(sort(tapply(comb.df$histological_type, comb.df$histological_type, length))))+
    ggtitle("GTEx data from SRAdb")
  ggsave(paste(outdir,"bodysite-bars.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(body_site,fill=sex))+geom_bar() + coord_flip() + theme_bw() + 
    scale_x_discrete(limits=names(sort(tapply(comb.df$body_site, comb.df$body_site, length))))+
    ggtitle("GTEx data from SRAdb")
  ggsave(paste(outdir,"histological_type-bars.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(bases))+ geom_density()+ theme_bw()+
    ggtitle("distribution of number of bases read in RNA-seq")
  ggsave(paste(outdir,"bases_per_run",sep=""),height=5,width=5)
  
  ggplot(comb.df, aes(x=body_site))+geom_bar()
  
  vals <- lapply(unique(comb.df$histological_type), function(x)dim(comb.df[which(comb.df$histological_type == x),])[1])
  names(vals) <- unique(comb.df$histological_type)
  names(vals[vals > 10])
  
  valsb <- lapply(unique(comb.df$body_site), function(x)dim(comb.df[which(comb.df$body_site == x),])[1])
  names(valsb) <- unique(comb.df$body_site)
  names(valsb[valsb > 20])
  valsb20 = valsb[valsb > 20]
  
  ggplot(comb.df[which(comb.df$body_site %in% names(valsb[valsb > 20])),], aes(x=body_site)) + 
    geom_bar() + coord_flip() + 
    thisTheme2 + ggtitle("GTEx RNA-seq\nbody sites with 20 or more samples") + 
    xlab("body site") + 
    scale_x_discrete(limits=names(valsb20)[order(as.numeric(valsb20))])
  
}
#reorder(Position,Position, function(x)-length(x)))

#library_selection
#library_strategy
#library_name


plotGTEx2014 <- function(){
  d <- getGTExAnnot()
  m <- readInGTExAllMeta()
  d$SAMPIDdash = gsub(d$SAMPID, pattern="\\.",replacement="-")
  m$gtexId= sapply(m$gtexId,toupper)
  d$SAMPIDdash= sapply(d$SAMPIDdash,toupper)
  comb <- merge(d,m, by.x="SAMPIDdash",by.y="gtexId")
  comb$subject_id <- unlist(sapply(comb$sample_attribute, function(x)do.call("[[",list(convertToList(x),"gap_subject_id"))),use.names=FALSE)
  sex.list <- sapply(comb$sample_attribute, function(x)do.call("[[",list(convertToList(x),"sex")))
  #body_site <- sapply(comb$sample_attribute, function(x)unlist(do.call("$",list(convertToList(x),"body_site")),
  #                                                             use.names=FALSE),simplify="vector")
  comb$body_site <- comb$SMTS
  comb$histological_type <- comb$SMTSD
  comb$gap_subject_id <- factor(comb$subject_id)
  
  indCount <- tapply(comb$gap_subject_id, comb$subject_id, length)
  body_siteCount <- tapply(comb$SMTSD, comb$SMTSD, length)
  
  
  ggplot(as.data.frame(indCount), aes(x=indCount))+geom_bar()+
    ggtitle("GTEX 2014\nNumber of samples per participant")
  ggsave(paste(outdir,"samplesPerPerson-heatmap-2014.pdf",sep=""),height=6,width=12)
  
  outdir <- getFullPath("plots/coexpr/GTEx")
  
  ggplot(comb, aes(as.numeric(gap_subject_id),SMTSD,fill=factor(SMTS)))+geom_tile() +theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb$SMTSD, comb$SMTSD, length))))+
    ggtitle("GTEx data from SRAdb\n2014")
  ggsave(paste(outdir,"bodysite-heatmap-2014.pdf",sep=""),height=6,width=12)
  
  
  
  bodySites30 <- names(table(comb$SMTSD)[which(table(comb$SMTSD) > 30)])
  subjCount10<- names(table(comb$gap_subject_id)[which(table(comb$gap_subject_id) > 9)])  
  comb30 <- comb[which(comb$SMTSD %in% bodySites30),]
  comb30$id <- as.numeric(factor(comb30$subject_id))
  
  ggplot(comb30, aes(x=as.numeric(factor(gap_subject_id)),y=SMTSD))+geom_tile() +theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb30$SMTSD, comb30$SMTSD, length))))+
    scale_x_discrete(limits=names(sort(tapply(comb30$gap_subject_id, comb30$gap_subject_id, length))),
                     labels=rep("",length.out=length(comb30$gap_subject_id)))+
    ggtitle("GTEx data from SRAdb\n2014\n30 or more samples present")+
    theme(axis.title.y = element_text(size = rel(1)))+
    theme(axis.title.x = element_text(size = rel(1)))+
    theme(axis.text.x = element_text(size = rel(0.6), angle = 90))+
    theme(axis.text.y = element_text(size = rel(0.6)))+
    xlab("Individual Donor") + ylab("Histological Site Sampled")
  ggsave(paste(outdir,"bodysite-heatmap-2014-above30examples.png",sep=""),height=4,width=7)
  
  brain <- names(table(comb$SMTSD)[grep(names(table(comb$SMTSD)),pattern="Brain")])
  testis <- names(table(comb$SMTSD)[grep(names(table(comb$SMTSD)),pattern="Testis")])
  
  comb30brainBalls <- comb[which(comb$SMTSD %in% unique(c(brain,bodySites30))),]
  ggplot(comb30brainBalls, aes(x=as.numeric(factor(gap_subject_id)),y=SMTSD))+geom_tile() +theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb30$SMTSD, comb30$SMTSD, length))))+
    scale_x_discrete(limits=names(sort(tapply(comb30$gap_subject_id, comb30$gap_subject_id, length))),
                     labels=rep("",length.out=length(comb30$gap_subject_id)))+
    ggtitle("GTEx data from SRAdb\n2014\n30 or more samples present")+
    theme(axis.title.y = element_text(size = rel(1)))+
    theme(axis.title.x = element_text(size = rel(1)))+
    theme(axis.text.x = element_text(size = rel(0.6), angle = 90))+
    theme(axis.text.y = element_text(size = rel(0.6)))+
    xlab("Individual Donor") + ylab("Histological Site Sampled")
  ggsave(paste(outdir,"bodysite-heatmap-2014-above30examples.png",sep=""),height=4,width=7)
  
  
  comb30s10 <- comb[which(comb$SMTSD %in% bodySites30 & comb$gap_subject_id %in% subjCount10),]
  ggplot(comb30, aes(gap_subject_id,SMTSD,fill=factor(SMTS)))+geom_tile() +theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb30$SMTSD, comb30$SMTSD, length))))+
    scale_x_discrete(limits=names(sort(tapply(comb30$gap_subject_id, comb30$gap_subject_id, length))))+
    ggtitle("GTEx data from SRAdb\n2014\n30 or more samples present")
  ggsave(paste(outdir,"bodysite-heatmap-2014-above30examplesIndividualFound10times.pdf",sep=""),height=6,width=12)
  
  
  cellTypes <-  combn(names(sort(tapply(comb30$SMTSD, comb30$SMTSD, length),decreasing=FALSE))[1:10],5)
  indWcellTypes <- apply(cellTypes,2,function(x)with(comb[which(comb$SMTSD %in% x),],sum(tapply(SMTSD,gap_subject_id,length) == 5,na.rm=TRUE)))
  cellTypes[,which(indWcellTypes == max(indWcellTypes))]
  
  
  
  pairsInd <-  combn(names(sort(tapply(comb$gap_subject_id, comb$gap_subject_id, length),decreasing=FALSE)),2)
  cellTypesInComman <- apply(pairsInd,2,function(x)with(comb[which(comb$gap_subject_id %in% x),],sum(tapply(gap_subject_id,SMTSD,length) == 2,na.rm=TRUE)))
  cellTypes[,which(indWcellTypes == max(indWcellTypes))]
  
  
  
  ggplot(comb.df, aes(histological_type,fill=sex))+geom_bar() + coord_flip() + theme_bw() +
    scale_x_discrete(limits=names(sort(tapply(comb.df$histological_type, comb.df$histological_type, length))))+
    ggtitle("GTEx data from SRAdb")
  ggsave(paste(outdir,"bodysite-bars.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(body_site,fill=sex))+geom_bar() + coord_flip() + theme_bw() + 
    scale_x_discrete(limits=names(sort(tapply(comb.df$body_site, comb.df$body_site, length))))+
    ggtitle("GTEx data from SRAdb")
  ggsave(paste(outdir,"histological_type-bars.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(bases))+ geom_density()+ theme_bw()+
    ggtitle("distribution of number of bases read in RNA-seq")
  ggsave(paste(outdir,"bases_per_run",sep=""),height=5,width=5)
  
  ggplot(comb.df, aes(x=body_site))+geom_bar()
  
  vals <- lapply(unique(comb.df$histological_type), function(x)dim(comb.df[which(comb.df$histological_type == x),])[1])
  names(vals) <- unique(comb.df$histological_type)
  names(vals[vals > 10])
  
  valsb <- lapply(unique(comb.df$body_site), function(x)dim(comb.df[which(comb.df$body_site == x),])[1])
  names(valsb) <- unique(comb.df$body_site)
  names(valsb[valsb > 20])
  valsb20 = valsb[valsb > 20]
  
  ggplot(comb.df[which(comb.df$body_site %in% names(valsb[valsb > 20])),], aes(x=body_site)) + 
    geom_bar() + coord_flip() + 
    thisTheme2 + ggtitle("GTEx RNA-seq\nbody sites with 20 or more samples") + 
    xlab("body site") + 
    scale_x_discrete(limits=names(valsb20)[order(as.numeric(valsb20))])
  
}


genWebsiteKey <- function(SRA="SRR613771"){
  system(paste("cd /data/wespisea/gtex/sraDB; ~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/test-sra",SRA,"2>&1 > ~/sandbox/sratoolTest &"))
  Sys.sleep(30)
}

getWebsiteKey <- function(){
  lines <- readLines(con="~/sandbox/sratoolTest")
  lines <- lines[grep(x=lines, "Remote http")]
  lines <- lines[grep(x=lines, "http://gap-upload.ncbi.nlm.nih.gov/")]
  
  
  #   
  # p <- pipe("cat ~/sandbox/sratoolTest |grep 'Remote http' |grep 'gap'")
  # res <- readLines(p)
  as.character(unlist(strsplit(x=as.character(unlist(strsplit(x=lines[1],split ="gap-upload.ncbi.nlm.nih.gov/")))[-1],split="/")[1]))[1]
}

generateFastqConvert <- function(vec){
  cmd<-paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump ",
              "-O /data/wespisea/gtex/fastq/ --split-files -gzip /data/wespisea/gtex/sra/",vec,".sra")
}


downloadFileMissing_url_getkey <- function(){
  df <- readInGTExAllMeta()
  df$fastqConvert <- generateFastqConvert(df$run_accession)
  
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$sraFile <- paste0("/data/wespisea/gtex/sra/",df$run_accession,".sra")
  df$haveFiles <- file.exists(df$fastq1)
  df$haveSra <- file.exists(df$sraFile)
  df.need <- df[which(df$haveSra== FALSE),]
  
  write(paste0(df.need$downloadCmd,rep(c(" &"," "),length=length(df.need$downloadCmd))), file="~/sandbox/downloadGTEx_need.sh") 
  # http://gap-upload.ncbi.nlm.nih.gov/E2004057-252A-4BC9-ADA3-2463E4AC9B98/SRR612551.sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47
  
  genWebsiteKey() #takes 30 seconds....
  key <- getWebsiteKey()
  print(paste("website key is ", key))
  o <- paste0("wget -O /data/wespisea/gtex/sra/",df.need$run_accession,".sra 'http://gap-upload.ncbi.nlm.nih.gov/",key,"/",
              df.need$run_accession,".sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47' --continue",
              rep(c(" &"," "),length=length(df.need$run_accession)))
  
  
  len <- length(o)
  midpoint <- floor(len/2)
  write(c(o[1:midpoint],df.need$fastqConvert[1:midpoint]), file = "~/sandbox/downloadGTEx2_wget1.sh")
  write(c(o[(midpoint+1):len],df.need$fastqConvert[(midpoint+1):len]), file = "~/sandbox/downloadGTEx2_wget2.sh")
  
  
  
  convert <- with(df,which(haveSra == TRUE & haveFiles == FALSE))
  df.convert <- df[convert,]
  o2 <- df.convert$fastqConvert
  
  len2 <- length(o2)
  midpoint2 <- floor(len2/2)
  write(o2[1:midpoint2], file = "~/sandbox/downloadGTEx2_convert1.sh")
  write(o2[(midpoint2+1):len2], file = "~/sandbox/downloadGTEx2_convert2.sh")
  
  
  s0 <- "find /data/wespisea/gtex/sra/ -name \"*.sra\" -size 0 -delete"
  s1 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_wget1.sh\""
  s2 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_wget2.sh\""
  s3 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_convert1.sh\""
  s4 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_convert2.sh\""
  cat(paste0(c(s0,s1,s2,s3,s4),collapse="\n"))
  c(s0,s1,s2,s3,s4)
  
}

transferToNearline <- function(localFile,remoteDir="/farline/umw_zhiping_weng/wespisea/gtex/sra/"){
  paste0("scp ", localFile, " aw30w@ghpcc06.umassrc.org:",remoteDir)
  
}

rmSraFoundInFarline <- function(sra,local="/data/wespisea/gtex/sra/",
                                remote="/farline/umw_zhiping_weng/wespisea/gtex/sra/"){
  sraLocal <- paste0(local,sra)
  sraRemote <- paste0(remote,sra)
  localExists <- file.exists(sraLocal)
  remoteExists <- sapply(sraRemote,hpc.file.exists)
  sraLocalRemote <- sraLocal[localExists & remoteExists]
  paste0("rm ",sraLocal)
}

rmSraFiles <- function(){
  df <- readInGTExAllMeta()
  df$fastqConvert <- generateFastqConvert(df$run_accession)
  
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$sraFile <- paste0("/data/wespisea/gtex/sra/",df$run_accession,".sra")
  df$sraFileBase <- paste0(df$run_accession,".sra")
  df$haveFiles <- file.exists(df$fastq1)
  df$haveSra <- file.exists(df$sraFile)
  rmCmd <- rmSraFoundInFarline(df$sraFileBase)
  
}

#run when all SRA's are downloaded & we need to transfer and convert
convertAndTransferSRA <- function(){
  df <- readInGTExAllMeta()
  df$fastqConvert <- generateFastqConvert(df$run_accession)
  
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$sraFile <- paste0("/data/wespisea/gtex/sra/",df$run_accession,".sra")
  df$sraFileBase <- paste0(df$run_accession,".sra")
  df$haveFiles <- file.exists(df$fastq1)
  df$haveSra <- file.exists(df$sraFile)
  
  df.need <- df[which(df$haveSra== FALSE),]
  
  #write(paste0(df.need$downloadCmd,rep(c(" &"," "),length=length(df.need$downloadCmd))), file="~/sandbox/downloadGTEx_need.sh") 
  # http://gap-upload.ncbi.nlm.nih.gov/E2004057-252A-4BC9-ADA3-2463E4AC9B98/SRR612551.sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47
  
  genWebsiteKey() #takes 30 seconds....
  key <- getWebsiteKey()
  print(paste("website key is ", key))
  o <- paste0("wget -O /data/wespisea/gtex/sra/",df.need$run_accession,".sra 'http://gap-upload.ncbi.nlm.nih.gov/",key,"/",
              
              df.need$run_accession,".sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47' --continue",
              rep(c(" &"," "),length=length(df.need$run_accession)))
  
  
  len <- length(o)
  midpoint <- floor(len/2)
  write(c(o[1:midpoint],o[1:midpoint]), file = "~/sandbox/downloadGTEx2_wget1.sh")
  write(c(o[(midpoint+1):len],o[(midpoint+1):len]), file = "~/sandbox/downloadGTEx2_wget2.sh")
  
  
  
  convert <- with(df,which(haveSra == TRUE & haveFiles == FALSE))
  df.convert <- df[convert,]
  o2 <- df.convert$fastqConvert
  
  len2 <- length(o2)
  midpoint2 <- floor(len2/2)
  write(o2[1:midpoint2], file = "~/sandbox/downloadGTEx2_convert1.sh")
  write(o2[(midpoint2+1):len2], file = "~/sandbox/downloadGTEx2_convert2.sh")
  
  transfer <- with(df,which(haveSra == TRUE & haveFiles == TRUE))
  df.transfer <- df[transfer,]
  o3 <- transferToNearline(df.transfer$sraFile)
  
  len3 <- length(o3)
  midpoint3 <- floor(len3/2)
  write(o3[1:midpoint3], file = "~/sandbox/downloadGTEx2_transfer1.sh")
  write(o3[(midpoint3+1):len3], file = "~/sandbox/downloadGTEx2_transfer2.sh")
  
  

  s0 <- "find /data/wespisea/gtex/sra/ -name \"*.sra\" -size 0 -delete"
  s1 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_wget1.sh\""
  s2 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_wget2.sh\""
  s3 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_convert1.sh\""
  s4 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_convert2.sh\""
  s5 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_transfer1.sh\""
  s6 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_transfer2.sh\""
  cat(paste0(c(s0,s1,s2,s3,s4,s5,s6),collapse="\n"))
  c(s0,s1,s2,s3,s4,s5,s6)
  
  d <- getGTExAnnot()
  m <- df
  d$SAMPIDdash = gsub(d$SAMPID, pattern="\\.",replacement="-")
  m$gtexId= sapply(m$gtexId,toupper)
  d$SAMPIDdash= sapply(d$SAMPIDdash,toupper)
  comb <- merge(d,m, by.x="SAMPIDdash",by.y="gtexId")
  exportAsTable(file="/data/wespisea/gtex/annot/GTExJan2014_wAnnot.tab",df=comb)
  
}

SRAstatusGood <- function(sraDir="/data/wespisea/gtex/sra/"){
  sraDir <- "/data/wespisea/gtex/sra/"
  # sraDir <- "/home/wespisea/sandbox/testSRA"
  p <- pipe(paste("ls",sraDir))
  files <- file.path(sraDir,readLines(p))
  close(p)
  sizeDistro <- sapply(files, function(x)file.info(x)$size)
  if (sum(sizeDistro == 0) > 10){
    FALSE
  }
  TRUE
}


runDownloadMoniter <- function(){
  
  while(TRUE){
    if (FALSE == SRAstatusGood()){
      cmds <- downloadFileMissing_url_getkey()
      system(cmds[1])
      Sys.sleep(10)
      system(cmds[2])
      system(cmds[3])
      Sys.sleep(10)
    } 
    Sys.sleep(3600)
  }
} 

runConvertToFastq <- function(){
  downloadFileMissing()
}
checkFastqFile <- function(){
  #/home/wespisea/bin/FastQC/fastqc -q -o ~/sandbox SRR1092397_1.fastq.gz
  # cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip /data/wespisea/gtex/sra/SRR1072247.sra &
}

