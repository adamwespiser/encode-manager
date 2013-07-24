library(ggplot2)


file <- "/Users/adam/Desktop/cellTfLab.tab"

df <- read.table(file=file,header=TRUE)
 df$one=1
 cdata<-cast(df, lab + cell ~ antibody,value="one")
df$lab == "HudsonAlpha" -> x
df$lab == "Stanford"    -> y
ha <- df[which(df$lab == "HudsonAlpha"),]
s <- df[which(df$lab == "Stanford"),]
df.1 <- rbind(ha,s)

ggplot(df.1,aes(x=cell,y=antibody))+geom_tile()+theme_bw()+labs(title="PeakSeq")+facet_wrap(~lab) +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
  axis.text.y = element_text(colour="grey20",size=7,angle=0,hjust=1,vjust=0,face="plain"))
 ggsave(file="~/Dropbox/wespiser_zlab_papers/TfPeakSeq/PeakSeq.pdf",height=14,width=6)

ggplot(df,aes(x=cell,y=antibody))+geom_tile()+theme_bw()+labs(title="PeakSeq")+facet_grid(.~lab) +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=7,angle=0,hjust=1,vjust=0,face="plain"))
ggsave(file="~/Dropbox/wespiser_zlab_papers/TfPeakSeq/PeakSeqAll.pdf",height=14,width=13)


ggplot(df,aes(x=cell,y=antibody))+geom_tile()+theme_bw()+labs(title="PeakSeq\nAll Combined") +
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=7,angle=0,hjust=1,vjust=0.5,face="plain"))
ggsave(file="~/Dropbox/wespiser_zlab_papers/TfPeakSeq/PeakSeqCombined.pdf",height=14,width=3)

ggplot(s,aes(x=cell,y=antibody))+geom_tile()+theme_bw()+labs(title="Stanford-PeakSeq")
ggsave(file="~/Dropbox/wespiser_zlab_papers/TfPeakSeq/StanfordTf-PeakSeq.pdb")

