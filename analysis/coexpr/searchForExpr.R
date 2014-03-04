library(SRAdb)
#sqlfile <- getSRAdbFile()

sqlfile <- "/home/wespisea/data/SRAmetadb.sqlite" 
sra_con <- dbConnect(SQLite(),sqlfile)
sra_tables <- dbListTables(sra_con)
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
  # get sample attributes where sample attrs has 16 members
  # sampleAttr.16   <- sapply(gtex.16.df$sample_attribute[1:100], function(x)as.data.frame(convertToList(x)))
  
  # convert sample attributes to df
  sampleAttr.df <- ldply(gtex.16.df$sample_attribute, function(x)as.data.frame(convertToList(x)))
  
  file.df <- ldply(gtex.16.df$run_accession, function(x)listSRAfile(x,sra_con))
  
  comb.df <- cbind(gtex.16.df,sampleAttr.df,file.df)
}

plotGTEx <- function(){
  comb.dfdf <- getGTExSra()
  outdir <- getFullPath("plots/coexpr/GTEx")
  
  ggplot(comb.df, aes(as.numeric(gap_subject_id),body_site,fill=sex))+geom_raster() +theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb.df$body_site, comb.df$body_site, length))))
  ggsave(paste(outdir,"bodysite-heatmap.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(as.numeric(gap_subject_id),histological_type,fill=sex))+geom_raster()+theme_bw()+ 
    scale_y_discrete(limits=names(sort(tapply(comb.df$histological_type, comb.df$histological_type, length))))
  ggsave(paste(outdir,"histological_type-heatmap.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(histological_type,fill=sex))+geom_bar() + coord_flip() + theme_bw() +
    scale_x_discrete(limits=names(sort(tapply(comb.df$histological_type, comb.df$histological_type, length))))
  ggsave(paste(outdir,"bodysite-bars.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(body_site,fill=sex))+geom_bar() + coord_flip() + theme_bw() + 
    scale_x_discrete(limits=names(sort(tapply(comb.df$body_site, comb.df$body_site, length))))
  ggsave(paste(outdir,"histological_type-bars.pdf",sep=""),height=12,width=5)
  
  ggplot(comb.df, aes(bases))+ geom_density()+ theme_bw()+
    ggtitle("distribution of number of bases read in RNA-seq")
  ggsave(paste(outdir,"bases_per_run",sep=""),height=12,width=5)
  
}
#reorder(Position,Position, function(x)-length(x)))

#library_selection
#library_strategy
#library_name
