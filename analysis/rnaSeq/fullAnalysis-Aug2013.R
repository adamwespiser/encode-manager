#The purpose of this script is to run: lncRNA comparison, PCA, logistic Regression, and eigenRank 
# the output will be directed to one directory
#
#
#

source("/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/eigenRank_08272013.R")
source("/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/logisticRegPcaExprData.R")

#LOG Reg 

runLogReg(outdir=outdir,lncPca=lncExpr.df)



#################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################



  ## EigenRank Setup
sourceFile = getFullPath("data/lncExprWithStats_transEachSample.tab")
df = readInTable(sourceFile)
df = df[which(df$averageExpr != 0),]
df[["gene_id_short"]] =  sapply(as.character(df$gene_id),function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})

doubleCols = colnames(df)[as.vector(sapply(colnames(df),function(x)(typeof(df[1,x]) == "double")))]
exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]

lncFound.df = getEnslist()
lncFound.df = unique(data.frame(ensembl_gene_id=lncFound.df$ensembl_gene_id,external_gene_id=lncFound.df$external_gene_id,gene_biotype=lncFound.df$gene_biotype))
lncFound.df = lncFound.df[which(lncFound.df$ensembl_gene_id %in% df[["gene_id_short"]]), ]

df$label = 0
df[which(df$gene_id_short %in% lncFound.df$ensembl_gene_id),"label"] = 1

biotype.df = readInTable(getFullPath("data/lnc_biotype.tab"))
df = merge(df,biotype.df[c("gene_id_short","bm_biotype")],by.x="gene_id_short",by.y="gene_id_short")


cols.list = list(lpa=exprCols.lpa,lnpa=exprCols.lnpa,bothPullDowns=c(exprCols.lpa,exprCols.lnpa))
biotypes.vec = c("antisense","lincRNA","processed_transcript","all_biotypes")
rows.list = list(antisense = df[which(df$bm_biotype == "antisense"),"gene_id_short"],
                   processed_transcript = df[which(df$bm_biotype == "processed_transcript"),"gene_id_short"],
                   lincRNA = df[which(df$bm_biotype == "lincRNA"),"gene_id_short"],
                   all_biotypes = df[c("gene_id_short")]
                   )

basedir = getFullPath("plots/rnaSeq-eigenRank/functionalTypes/")


## LogReg Setup
lnc.pca.df <- getPcaData()

#put label y=1 at top...
lnc.pca.df <- lnc.pca.df[order(lnc.pca.df$label,decreasing=TRUE),]




for(columns in c("lpa","lnpa","bothPullDowns")){
  expr.cols = cols.list[[columns]]
  for(biotype in biotypes.vec){
   
    ## EigenRank
    expr.rows = rows.list[[biotype]]
    print(paste("starting",columns,biotype))
    outdir = paste(basedir,columns,"-",biotype,sep="")
    print(outdir)
    if(!file.exists(outdir)){dir.create(outdir)}
    #dir.create(outdir)
    local.df = df[which(df$gene_id_short %in% as.vector(expr.rows)),]
    nolab = local.df[which(local.df$label == 0),]
    lab = local.df[which(local.df$label == 1),]
    outdir = paste(outdir,"/",sep="")
    
    print("egien")
    plotEigenVectorsDensity(lab=lab,nolab=nolab,outdir=outdir,
                            filename=paste(columns,biotype,sep="-"),cols=expr.cols ,titleMsg="")
    
    
    ## Log. Regression
    
    
  }
}




#function(df,lab,nolab,outdir,filename,cols,titleMsg="")





