# Purpose: compare the features of of known lncRNA and the rest of lncRNA
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/createMLdataset.R")

library(Rcpp)
library(inline)
library(RcppArmadillo)



home <- Sys.getenv("HOME")

projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

#boiler plate helpers
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)

#pull in exprLib for some helper functions, and getENSGfromBiomartByRefseq.R for the gene BioMart grabs..
#source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R"#)
#source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R")

srcCalcEigenCpp <- '
  Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericVector yr(ys);
int n = Xr.nrow(),  k = Xr.ncol();
arma::mat X(Xr.begin(), n, k, false);
arma::colvec y(yr.begin(), yr.size(), false);
arma::colvec y_prev = y;
arma::colvec y_temp = y;


arma::mat sim =  (X*trans(X));
sim = arma::pow(sim,-1);

int i = 0;
double diff = 10000000.0;
double y_norm = 0;
while(diff > 1e-10 && i < 2000){
i++; 
y_prev = y;

y_temp = sim * y;
y_norm = norm(y_temp,2);
y = y_temp / y_norm ;
//diff = sum(abs(y_prev - y)) ;
diff = norm(y_prev - y, 2);
} 

return Rcpp::List::create(Rcpp::Named("y") = y,
Rcpp::Named("ytemp") = y_temp,
Rcpp::Named("yprev") = y_prev,
Rcpp::Named("ynorm") = y_norm,
Rcpp::Named("converge")    = abs(sum(y_prev - y)),
Rcpp::Named("iters")       = i,
Rcpp::Named("diff")       = diff);

'  # end of src4....
# ... 
CalcEigenCpp = cxxfunction(signature(Xs="numeric", ys="numeric"), srcCalcEigenCpp, plugin="RcppArmadillo")

# all three of these methods give the same rank...
srcCalcEigenCppD <- '
Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericVector yr(ys);
int n = Xr.nrow(),  k = Xr.ncol();
arma::mat X(Xr.begin(), n, k, false);
arma::colvec y(yr.begin(), yr.size(), false);
arma::colvec y_prev = y;
arma::colvec y_temp = y;

arma::mat sim =  X*trans(X);
arma::mat distance = sim;
double found_zero = 0;



int i;

for(i = 0; i < n; i++){
double bii = sim(i,i);
int j;
for (j = i ; j < n; j++){
if(i == j){
distance(i,j) = 0;
continue;
} 
double distValue = 0;
double bjj = sim(j,j);
double bij = sim(i,j);
double dist2 = bii + bjj - 2*bij; 
if (pow(dist2,2) < 1e-15){
distValue = 0;
}  
else{
distValue = sqrt(dist2);

}
double simValue = 1 / distValue;
if (std::isnan(simValue) || std::isinf(simValue)){
simValue = 0;

}
distance(i,j) = distance(j,i) = simValue;

}
}

i = 0;
double diff = 10000000.0;
double y_norm = 0;
while(diff > 1e-15 && i < 2000){
i++; 
y_prev = y;

y_temp = distance * y;
y_norm = norm(y_temp,2);
y = y_temp / y_norm ;
diff = norm(y_prev - y, 2);
} 

arma::colvec y_norm_out = (y / sum(y)) * n;

return Rcpp::List::create(Rcpp::Named("y") = y,
Rcpp::Named("ytemp") = y_temp,
Rcpp::Named("yprev") = y_prev,
Rcpp::Named("yn") = y_norm_out,
Rcpp::Named("ynorm") = y_norm,
Rcpp::Named("n")    =n, 
Rcpp::Named("iters")       = i,
Rcpp::Named("foundZero")       = found_zero,
Rcpp::Named("diff")       = diff,
Rcpp::Named("distM")  = distance,
Rcpp::Named("sim")   = sim);

' 
# end of src4....
# ... 

CalcEigenCppD = cxxfunction(signature(Xs="numeric", ys="numeric"), srcCalcEigenCppD, plugin="RcppArmadillo")
# CalcEigenCpp(Xs=(m+ 1),y=ystart) -> a1
# CalcEigenCpp(Xs=(m+ 100),y=ystart) -> a2
# CalcEigenCpp(Xs=(m+ 1000),y=ystart) -> a3




plotEigenVectorsPvals = function(nolab,lab,outdir,filename,cols,titleMsg=""){
 
  if(!file.exists(outdir)){dir.create(path=outdir,recursive=TRUE)}
  eigenRankOutfile = paste(outdir,"/",filename,"rank.pdf",sep="")
  pvalueOutfile = paste(outdir,"/",filename,"pvalues.pdf",sep="")
  fullEigenOutfile= paste(outdir,"/",filename,"rank_all.pdf",sep="")
  rankVexprOutfile= paste(outdir,"/",filename,"rank_all-vs-aveExpr.pdf",sep="")
  title = paste("PageRank of lncRNA expr\nfunctional vs rest of lncRNA",titleMsg,sep="\n") 
  
  
  eigen.df = ddply(nolab, .(group), function(x)calcEigenRankLabelNoLabel(lab = lab,nolab = x,cols = cols))
  eigenpvals.df = ddply(eigen.df, .(group), function(x)data.frame(pval = wilcox.test(x[x$type == 0,"rank"], x[x$type == 1,"rank"])[["p.value"]]))
  
  
  # colnames(eigen.df) = c("group", "rank", "type", "lncRnaName")
  ggplot(eigen.df, aes(x = log(rank), fill = factor(type))) + geom_density(alpha=I(0.6)) + facet_wrap(~group)+theme_bw()+ 
    ggtitle(title)
  ggsave(file=eigenRankOutfile)
  
  ggplot(eigenpvals.df,aes(x=group,y=log10(pval)))+geom_bar(binwidth=1,stat="identity")+coord_flip()+
    ggtitle(title)
  ggsave(file=pvalueOutfile)
  
  
  #full eigenVector calc:
  
  # exprCols.lnpa # exprCols.lpa # maxCols
  mat.df = rbind(nolab[,cols],lab[,cols])
  mat = (as.matrix(mat.df)) + 1
  initialGuess = rep(1,dim(mat)[1])
  e.out = CalcEigenCpp(Xs = (mat + 1), y = runif(dim(mat.df)[1]))
  mat.df$rank = e.out$y
  mat.df$label = c(nolab$label,lab$label)
  mat.df$lncRnaName = c(nolab$lncRnaName,lab$lncRnaName)
  mat.df$rank = (mat.df$rank / sum(mat.df$rank))*length(mat.df$rank)
  mat.df$averageExpr = c(nolab$averageExpr, lab$averageExpr)
  m <<- mat.df
  
  ggplot(mat.df,aes(x=rank,fill=factor(label)))+geom_density(alpha=I(0.6))+theme_bw() + xlim(0.8,1)+
    ggtitle(title)
  ggsave(file=fullEigenOutfile)
  
  ggplot(mat.df,aes(y=averageExpr,x=rank,color=factor(label)))+geom_point()+theme_bw() +
    ggtitle("Rank versus average Expression")
  ggsave(file=rankVexprOutfile)
  
}


plotEigenVectorsDensity = function(lab,nolab,outdir,filename,cols,titleMsg=""){
  
  countNoLabel = dim(nolab)[1]
  countLabel = dim(lab)[1]
    
  if(!file.exists(outdir)){dir.create(path=outdir,recursive=TRUE)}
  eigenRankOutfile = paste(outdir,"/",filename,"rank.pdf",sep="")
  pvalueOutfile = paste(outdir,"/",filename,"pvalues.pdf",sep="")
  fullEigenOutfile= paste(outdir,"/",filename,"rank_all.pdf",sep="")
  fullEigenOutfile.log= paste(outdir,"/",filename,"rank_log.pdf",sep="")
  fullEigenOutfile.log2= paste(outdir,"/",filename,"rank_log-bars.pdf",sep="")
  fullEigenOutfile.log3= paste(outdir,"/",filename,"rank-bars.pdf",sep="")
  fullEigenOutfile.log4 = paste(outdir,"/",filename,"rank_log_abs.pdf",sep="")
  
  rankVexprOutfile= paste(outdir,"/",filename,"rank_all-vs-aveExpr.pdf",sep="")
  title = paste("PageRank of lncRNA expr\nfunctional vs rest of lncRNA",titleMsg,"\n",
                "count label=",countLabel," count nolabel=",countNoLabel, "\n",titleMsg,sep="") 
  
  #full eigenVector calc:
  # exprCols.lnpa # exprCols.lpa # maxCols
  mat.df = rbind(nolab[,cols],lab[,cols])
  mat = (as.matrix(mat.df)) + 1
  initialGuess = rep(1,dim(mat)[1])
  e.out = CalcEigenCppD(Xs = (mat + 1), y = runif(dim(mat.df)[1]))
  #e.out = CalcEigenCppD(Xs = (mat ), y = runif(dim(mat.df)[1]))
  #mat.df$rank = e.out$y_norm
  mat.df$rank = e.out$yn
  mat.df$label = c(nolab$label,lab$label)
  mat.df$lncRnaName = c(nolab$lncRnaName,lab$lncRnaName)
  mat.df$rank = (mat.df$rank / sum(mat.df$rank))*length(mat.df$rank)
  mat.df$averageExpr = c(nolab$averageExpr, lab$averageExpr)
  
  mat.df$logRank = log(mat.df$rank + 1)
  rank.sd =apply(mat.df$logRank, 2, sd)
  rank.median = median(mat.df$logRank)
  timesSd = 1
  xlim.range = c((rank.median - timesSd * rank.sd),(rank.median +  timesSd *rank.sd))
  plot.points = which(mat.df$logRank < xlim.range[2] & mat.df$logRank > xlim.range[1])
  
  rankn.sd =apply(mat.df$rank, 2, sd)
  rankn.median = median(mat.df$rank)
  timesSd = 1
  xlimn.range = c((rankn.median - timesSd * rankn.sd),(rankn.median +  timesSd *rankn.sd))
  plotn.points = which(mat.df$rank < xlimn.range[2] & mat.df$rank > xlimn.range[1])
  
  ggplot(mat.df[plot.points,],aes(x=logRank,fill=factor(label)))+geom_density(alpha=I(0.6))+theme_bw() + 
    xlim(xlim.range) +
    ggtitle(paste(title,"median +/- 2 ds",sep="\n") )
  ggsave(file=fullEigenOutfile.log)
  
  ggplot(mat.df ,aes(x=log(abs(rank)),fill=factor(label)))+geom_density(alpha=I(0.6))+theme_bw() + 
    ggtitle(paste(title,"all points",sep="\n") )
  ggsave(file=fullEigenOutfile.log4)
  
  ggplot(mat.df[plot.points,],aes(x=logRank,fill=factor(label))) + 
    geom_bar(binwidth=(xlim.range[2] - xlim.range[1])/30)+ 
    facet_wrap(~label,ncol=1,scale="free_y")+
    theme_bw() + 
    xlim(xlim.range) +
    ggtitle(paste(title,"median +/- 2 ds",sep="\n") )
  ggsave(file=fullEigenOutfile.log2)
  
  ggplot(mat.df[plotn.points,],aes(x=rank,fill=factor(label))) + 
    geom_bar(binwidth=((xlimn.range[2]) - (xlimn.range[1]))/30)+ 
    # geom_bar()+    
    facet_wrap(~label,ncol=1,scale="free_y")+
    theme_bw() + 
   xlim(xlimn.range) +
    ggtitle(paste(title,"median +/- 2 ds",sep="\n") )
  ggsave(file=fullEigenOutfile.log3)
  
  ggplot(mat.df,aes(x=rank,fill=factor(label)))+geom_density(alpha=I(0.6))+theme_bw() + 
    xlim(-0.1,1.1) +
    ggtitle(title)
  ggsave(file=fullEigenOutfile)
  
  ggplot(mat.df,aes(y=averageExpr,x=rank,color=factor(label)))+geom_point()+theme_bw() +
    ggtitle("Rank versus average Expression")
  ggsave(file=rankVexprOutfile)
  
}



eigenRankSplitUpIteration <- function(folder = getFullPath("plots/rnaSeq-eigenRank"), rankByCol = NULL,
                                      sourceFile = getFullPath("data/lncRnaExpr_ml.tab") ){
  df = readInTable(sourceFile)
  df = df[which(df$averageExpr != 0),]

  #divide dataset into labelled and unlabeled entries
  label.df   = df[which(df$label == 1),] 
  nolabel.df = df[which(df$label == 0),] 
  
  # (optional) sort by a column of coice
  nolabel.df = nolabel.df[order(nolabel.df[["averageExpr"]]), ]
  nolabel.df$index = seq_along(nolabel.df$label)
  
  
  #if (!is.null(rankByCol) == TRUE && (rankByCol %in% colnames(nolabel.df)) ){
  #  print(paste("using ", rankByCol, "column"))
  #  nolabel.df = nolabel.df[order(nolabel.df[[rankByCol]]), ]
#    nolabel.df$index = seq_along(nolabel.df$label)
#  }
  
  nolabel.df$group = cut(nolabel.df$index,breaks=pretty(nolabel.df$index,12))
  doubleCols = intersect(colnames(nolabel.df),
                         colnames(label.df))[as.vector(sapply(intersect(colnames(nolabel.df),colnames(label.df)),
                                                              function(x)(typeof(nolabel.df[1,x]) == "double")))]
  
  exprCols.lnpa = doubleCols[grep("longNonPolyA$",doubleCols)]
  exprCols.lpa  = doubleCols[grep("longPolyA$",doubleCols)]
  
  topColsToTake = 10
  rankedColsByPvalue = order(sapply(doubleCols,function(x)wilcox.test(nolabel.df[,x],label.df[, x])[["p.value"]]))
  rankedByDiff = order(sapply(doubleCols,
                              function(x){mean((nolabel.df[,x]- mean(df[,x]))/sd(df[,x])) 
                                          - mean((label.df[,x]- mean(df[,x]))/sd(df[,x]))}))
  
  maxCols =(doubleCols[rankedByDiff])[1:topColsToTake]
  
  
  
 
  plotEigenVectorsPvals(nolab = nolabel.df, lab = label.df, outdir = folder ,
                        filename = "eigenRank-rnaExpr-allcols-",
                        cols = which(colnames(nolabel.df) %in% c(exprCols.lnpa,exprCols.lpa)),
                        titleMsg = "all RNA expr experiments used")
  
  plotEigenVectorsPvals(nolab = nolabel.df, lab = label.df, outdir = folder,
                        filename = "eigenRank-longNonPolyA-allcols-",
                        cols = which(colnames(nolabel.df) %in% exprCols.lnpa),
                        titleMsg = "longNonPolyA RNA expr experiments used")
  
  plotEigenVectorsPvals(nolab = nolabel.df,lab = label.df, outdir = folder,
                        filename = "eigenRank-longPolyA-allcols-",
                        cols = which(colnames(nolabel.df) %in% exprCols.lpa),
                        titleMsg = "longPolyA RNA expr experiments used")
  

  plotEigenVectorsPvals(nolab = nolabel.df,lab = label.df,outdir = folder,
                        filename = "eigenRank-top10cols-allcols-",
                        cols = which(colnames(nolabel.df) %in% maxCols),
                        titleMsg = "top10 most different cols used")
}


main = function(){
  eigenRankSplitUpIteration(folder = getFullPath("plots/rnaSeq-eigenRank/randomGroups/"), rankByCol = NULL)
  eigenRankSplitUpIteration(folder = getFullPath("plots/rnaSeq-eigenRank/aveExprGroups/"), rankByCol = "averageExpr")
  
}



