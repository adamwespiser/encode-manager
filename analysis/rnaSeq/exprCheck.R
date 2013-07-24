df <- read.csv(file="~/Desktop/gencodeV7lncRNAexpr.csv")
lnc.df <- read.csv(file="~/Desktop/Gencode_lncRNAsv7_summaryTable_05_02_2012.csv")
lnc1.df <- read.table(file="/Users/adam/Desktop/Gencode_lncRNAsv7_summaryTable_05_02_2012_lncExression.tab")
library(plyr)
library(ddply)
library(reshape)
oneZero <- function(...){x = c(...)[1];y=c(...)[2];ifelse( (x*y == 0) && (x + y > 0),TRUE,FALSE )}

old <- df[c("filename","HELAS3")]
new <- lnc1.df[c("LncRNA_Txid","HELAS3")]
comb<-merge(old,new,by.x="filename",by.y="LncRNA_Txid")
colnames(comb) <- c("transcript_id", "derrien","rnaExprHub")
comb[which(apply(comb[2:3],1,oneZero)),] -> HELAS3




old <- df[c("filename","GM12878")]
new <- lnc1.df[c("LncRNA_Txid","GM12878")]
comb<-merge(old,new,by.x="filename",by.y="LncRNA_Txid")
colnames(comb) <- c("transcript_id", "derrien","rnaExprHub")
comb[which(apply(comb[2:3],1,oneZero)),] -> GM12878


old <- df[c("filename","HEPG2")]
new <- lnc1.df[c("LncRNA_Txid","HEPG2")]
comb<-merge(old,new,by.x="filename",by.y="LncRNA_Txid")
colnames(comb) <- c("transcript_id", "derrien","rnaExprHub")
comb[which(apply(comb[2:3],1,oneZero)),] -> HEPG2

old <- df[c("filename","K562")]
new <- lnc1.df[c("LncRNA_Txid","K562")]
comb<-merge(old,new,by.x="filename",by.y="LncRNA_Txid")
colnames(comb) <- c("transcript_id", "derrien","rnaExprHub")
comb[which(apply(comb[2:3],1,oneZero)),] -> K562

old <- df[c("filename","H1HESC")]
new <- lnc1.df[c("LncRNA_Txid","H1HESC")]
comb<-merge(old,new,by.x="filename",by.y="LncRNA_Txid")
colnames(comb) <- c("transcript_id", "derrien","rnaExprHub")

comb[which(apply(comb[2:3],1,oneZero)),] -> H1HESC
ggsave(file="~/Desktop/temp.pdf",width=7,height=7)
