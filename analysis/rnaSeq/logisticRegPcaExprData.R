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
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R")
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R")
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/exprLib.R")
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R")
source( "/Users/adam/work/research/researchProjects/encode/encode-manager/analysis/RnaSeq/pcaAnalysisWOoutlier.R")


getPcaData = function(exprCols=2:33){
  combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  combined.in.file = getFullPath("data/combinedExprWithStats_transEachSample.tab")
  lnc.in.file = getFullPath("data/lncExprWithStats_transEachSample.tab")
  
  lnc.colIndex = exprCols
  lnc.expr=readInTable(lnc.in.file)
  lnc.expr$gene_id_short <- sapply(lnc.expr$gene_id,function(x){as.vector(strsplit(x,"\\.")[[1]])[1]})
  func.df <- getEnslist()
 
  
  normalizePcaFactorsFromLoadingForRatioTestLocal <- function(pca,lncDataFrame,exprCols){
    lncDataFrame[exprCols] = transformDataByPCA(df=lncDataFrame[,exprCols],pca=pca)
    colnames(lncDataFrame)[exprCols] = paste0("Comp.",seq_along(exprCols))
    lncDataFrame
  }
  
  
  func.df <- within( func.df, {
    gene_id_short = ensembl_gene_id
    gene_id = gene_id_short
    lncRnaName = external_gene_id})
  f.df <- getAnnotLncDf(lncDf=lnc.expr,annotDf=func.df,exprCol=2:33,annotColName="lncRnaName")
  f.df[["lncRnaName"]] = ifelse(f.df$lncRnaName == "notFound",f.df$gene_id_short,f.df$lncRnaName)
  
  r1.vec = c("MALAT1","H19","RP11-255B23.3")
  r2.vec = c("ENSG00000235162","ENSG00000228474","ENSG00000175061","DANCR","SNGG1","ZFAS1","ENSG00000249790",
             "ENSG00000256329","ENSG00000235162","ENSG00000228474","ENSG00000249532","ENSG00000249502","MALAT1","H19","RP11-255B23.3")
  
  
  comb.df <- getAnnotLncDf(lncDf=lnc.expr,annotDf=f.df,exprCol=2:33,annotColName="lncRnaName")
  lncDf.reduced.2 <- f.df[which(!f.df$lncRnaName %in% r2.vec),]
  
  lnc.pca             <- princomp(lncDf.reduced.2[exprCols],cor=FALSE)

  exprCols=2:33
  lnc.pca.factors.r2.full <- normalizePcaFactorsFromLoadingForRatioTestLocal(pca=lnc.pca,lncDataFrame=comb.df,exprCols=exprCols)
  
  
  
  pca.df <- lnc.pca.factors.r2.full[exprCols]
  pca.df[["label"]] = ifelse(f.df[["withinSubset"]] == "false", 0 , 1) 
  pca.df
}

sigmoid <- function(x){ 1 / (1 + exp(-1 * x))}

#must accept theata as vector
costFunction <- function(X,y,thetaVec,k){
  theta = matrix(thetaVec,k,k)
  xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
  sig = sigmoid(xTx)
  #z = -y * log(sig) - (1 - y)*log(1 - sig)
  z1 = ifelse(y == 1, -log(sig),-log(1-sig))
  z2 = ifelse(z1 == Inf, 0, z1)
}

cfLambda <- function(X,y,k){
  function(thetaVec){
  theta = matrix(thetaVec,k)
  xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
  sig = abs(sigmoid(xTx))
  z1 = ifelse(y ==1, -log(sig),-log(1- sig))
  -1/length(y) * sum(ifelse(z1 == Inf, 0, z1))
  }
}

cfLambdaIterate <- function(X,y,k){
    function(thetaVec){
      theta = matrix(thetaVec,k)
      xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
      sig = sigmoid(xTx)
      cost = 0
      for( row in 1:length(y)){
        tmp = y[row]*log(sig[row]) + (1 - y[row])*log(1- sig[row])
        
      }
    cost * (-1/length(y))
}}


gradFunction <- function(X,y,thetaVec){
  theta = matrix(thetaVec,k)
  xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
  sig = sigmoid(xTx)
  s1 = sum(y - sig)
  (1*length(y))*as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) + s1
}

gfLambda <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
    sig = sigmoid(xTx)
   # print(length(y),length(sig))
    s1 = sum(sig-y)
    (1/length(y))* (as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1)
  }
  
}

gfLambdaIterate <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
    sig = sigmoid(xTx)
    outMat = matrix(rep(0,k*k),k)
    for( row in 1:length(y)){
      xRow = X[row,]
      rowSum = y[row] - sig[row]
      xxT =  X[row,] %*% t(X[row,]) * rowSum * 1/length(y)      
      for(j in 1:k){
       for(i in 1:k){ 
        outMat[i,j] = outMat[i,j] + xxT[i,j]
  
       }
       }
    }
    
    outMat 
    
  }
  
}


filterXbyOptimTheta <- function(X,thetaVec,k){
  theta = matrix(thetaVec,k)
  guess = apply(X,1,function(row) t(row) %*% matrix(o$par,k) %*% row)
  
  optim(fn=cfLamda(X,y,k),par=runif(k*k),method= "CG") -> o #$;o$value
   X = as.matrix(lnc.pca.df[1:n,1:k-1])
   X = cbind(rep(1,dim(X)[1]),X)
  y = as.matrix(lnc.pca.df[1:n,33])
  
  guess = apply(X,1,function(row) t(row) %*% matrix(o$par,k) %*% row)
  data.frame(x=guess,ellipse=guess,y=y,label=y) -> df
  ggplot(df,aes(x=log(x),fill=factor(y)))+geom_density(alpha=I(0.6))
  
  plotDistanceAwayRatio(df,label=y) -> df.plot
  ggplot(df.plot,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
}


filterXbyOptimTheta_getGuess <- function(X,thetaVec,k){
  theta = matrix(thetaVec,k)
  guess = apply(X,1,function(row) t(row) %*% matrix(o$par,k) %*% row)
  
 # optim(fn=cfLambda(X,y,k),par=runif(k*k),method= "CG") -> o #$;o$value
  X = as.matrix(lnc.pca.df[1:n,1:k-1])
  X = cbind(rep(1,dim(X)[1]),X)
  y = as.matrix(lnc.pca.df[1:n,33])
  
  guess = apply(X,1,function(row) t(row) %*% matrix(o$par,k) %*% row)
  guess
}


plotDistanceAwayRatio <- function(df,label,exprCols,
                                  df.center= suppressMessages(melt(ddply(df[exprCols],.(),numcolwise(mean)))[["value"]]),rand=FALSE,ellipse=TRUE){
  
  df[["withinSubset"]] = label
  if(rand == TRUE){
    df[["rand"]] = runif(length(df[["withinSubset"]]))
    df[["withinSubset"]] = df[order(df[["rand"]]),"withinSubset"]
  }
  
  if(ellipse == TRUE){
    df[["dist"]] = df[["ellipse"]]
    
  }
  else{
    df[["dist"]] = apply(df[exprCols],1,function(x) sqrt(sum((x - df.center)^2)))
    
  }
  dfo = df[order(max(df$dist) - df$dist),] # order from highest to lowest
  dfo[["TP"]] = cumsum(dfo$withinSubset)
  dfo[["predict"]] = seq_along(dfo$withinSubset)
  dfo[["P"]] = sum(dfo$withinSubset)
  dfo[["totalMembers"]] = length(dfo$withinSubset)
  dfo[["N"]] = dfo[["totalMembers"]] - dfo[["P"]]
  dfo.stats =within(dfo,{
    FP = predict - TP
    TN = (totalMembers - predict) - (P - TP)
    FN = P - TP
    sens = TP / (TP + FN)
    prec = TP / (TP + FP)
    TPR = TP / P
    FPR = FP / N
    FNR = FN / P
    TNR = TN / N
    Fscore = ( 2 * prec * sens ) / ( prec + sens )
    accuracy = (TP + TN) / (P + N)
    errorRate = (FP + FN) / (P + N)
    label = withinSubset
  })
  dfo.stats[c("predict","totalMembers","P","N","TP","FP","TN","FN","sens","prec","TPR","FPR","FNR","TNR","Fscore","accuracy","errorRate","dist","label")]
}

main <- function(){
  # load in data
  lnc.pca.df <- getPcaData()
  
  #put label y=1 at top...
  lnc.pca.df <- lnc.pca.df[order(lnc.pca.df$label,decreasing=TRUE),]
  
  # k is the number of features we are going to take
  # n is the number of lncRNA examples we are going to take, note that label=1 is sorted to the top
  k = 8 
  n = 400
  yCount = length(which(lnc.pca.df$label == 1))
  dataSample = c(1:94,sample(94:length(lnc.pca.df$label),(n-yCount)))
  dataSample = 1:n
  X = as.matrix(lnc.pca.df[dataSample,1:k-1])
  X = cbind(rep(1,dim(X)[1]),X)
  y = as.matrix(lnc.pca.df[dataSample,33])
    
  X.full = as.matrix(lnc.pca.df[,1:k-1])
  X.full = cbind(rep(1,dim(X.full)[1]),X.full)
  y.full = as.matrix(lnc.pca.df[,33])
  
  
  
  thetaGuess = runif(k*k)
  zeroGuess  = rep(0,k*k)
    
  #opt.lnc = optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=thetaGuess,method= "L-BFGS-B")
  #optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=thetaGuess,method= "CG") -> o #400
  #optim(fn=cfLamda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS")
  o <-optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS")
  # create cost function
  # create gradient
  # run optimizer
  
  results.df <- data.frame(trial=1:100,cost=rep(0,100))
  best.cost = 1000000
  best.theta =matrix(runif(k*k),k)
  for(i in 1:20){
    print(i)
    local =  optim(fn=cfLamda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS") 
    results.df$trial[i] = i
    results.df$cost[i] = local$value
    if (local$value < best.cost){
      best.cost = local$value
      best.theta = local$par
    }
    }
  
  dfo = plotDistanceAwayRatio(df=as.data.frame(X.full),label=y.full,exprCols=1:k,ellipse=FALSE)
  #ggplot(dfo,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  dfo$type = "circle"
  
  guess = apply(X.full,1,function(row) t(row) %*% matrix(best.theta$par,k) %*% row)
  df = data.frame(x=guess,ellipse=guess,y=y.full,label=y.full)
 # ggplot(df,aes(x=log(x),fill=factor(y)))+geom_density(alpha=I(0.6))
  df.plot= plotDistanceAwayRatio(df,label=y.full)
  df.plot$type = "ellipse"
  
  sample(y.full,length(y.full)) 
  df.random= plotDistanceAwayRatio(df,label=sample(y.full,length(y.full)))
  df.random$type = "randomLabel"
  
  #ggplot(df.plot,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  
  df.comb = rbind(dfo,df.plot)
  ggplot(df.comb,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()+facet_wrap(~type)
  ggsave("~/Desktop/ellipse1.pdf",height=4,width=6)
  
  df.comb.roc = rbind(dfo,df.plot,df.random)
  ggplot(df.comb.roc,aes(x=FPR,y=TPR,color=type))+geom_line()+geom_point()+theme_bw()+geom_abline(slope=1)+
    ggtitle("ROC Curve of radius based prediction on pca\n")+xlab("False Positive Rate")+ylab("True Positive Rate")
  ggsave(file="~/Desktop/ellipse-roc1.pdf",height=4,width=5)
  
  
}

runLogReg = function(lncDf,outdir,cols){
  k = length(cols) + 1
  n = dim(lncDf)[1]
  yCount = length(which(lncDf$label == 1))
  #dataSample = c(1:94,sample(94:length(lnc.pca.df$label),(n-yCount)))
 # dataSample = 1:n
  X = as.matrix(lncDf[,cols])
  X = cbind(rep(1,dim(X)[1]),X)
  y = as.matrix(lncDf[,"label"])
  
  
  
  
  thetaGuess = runif(k*k)
  zeroGuess  = rep(0,k*k)
  
  #opt.lnc = optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=thetaGuess,method= "L-BFGS-B")
  #optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=thetaGuess,method= "CG") -> o #400
  #optim(fn=cfLamda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS")
  o <-optim(fn=cfLambda(X,y,k),gr=gfLambda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS")
  # create cost function
  # create gradient
  # run optimizer
  
  results.df <- data.frame(trial=1:100,cost=rep(0,100))
  best.cost = 1000000
  best.theta =matrix(runif(k*k),k)
  for(i in 1:20){
    print(i)
    local =  optim(fn=cfLambda(X[sample(seq_along(X[,1]),50),],y,k),par=matrix(runif(k*k),k),method= "CG") 
    results.df$trial[i] = i
    results.df$cost[i] = local$value
    if (local$value < best.cost){
      best.cost = local$value
      best.theta = local$par
    }
  }
 # guess = filterXbyOptimTheta_getGuess(X,best.theta,k)
  
  
  
  dfo.circle = plotDistanceAwayRatio(df=as.data.frame(X),label=y,exprCols=cols,ellipse=FALSE)
  #ggplot(dfo,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  dfo.circle$type = "circle"
  
  guess = filterXbyOptimTheta_getGuess(X,best.theta,k)
  df.ellipse.tmp = data.frame(x=guess,ellipse=guess,y=y,label=y)
  # ggplot(df,aes(x=log(x),fill=factor(y)))+geom_density(alpha=I(0.6))
  df.ellipse = plotDistanceAwayRatio(df.ellipse.tmp,label=y,exprCols=cols,ellipse=TRUE)
  df.ellipse$type = "ellipse"
  
  #sample(y.full,length(y.full)) 
  #df.random= plotDistanceAwayRatio(df,label=sample(y.full,length(y.full)))
  #df.random$type = "randomLabel"
  
  #ggplot(df.plot,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  
  df.comb = rbind(dfo.circle,df.ellipse)
  ggplot(df.comb,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()+facet_wrap(~type)
  ggsave("~/Desktop/ellipse.pdf")
  
  ggsave(paste(outdir,"elipse.pdf"),height=4,width=6)
  
  df.comb.roc = rbind(dfo,df.plot,df.random)
  ggplot(df.comb.roc,aes(x=FPR,y=TPR,color=type))+geom_line()+geom_point()+theme_bw()+geom_abline(slope=1)+
    ggtitle("ROC Curve of radius based prediction on pca\n")+xlab("False Positive Rate")+ylab("True Positive Rate")
  ggsave(paste(outdir,"elipse-ROC.pdf"),height=4,width=5)
  
  
}





mainTest <- function(){
  k = 3 
  n = 100
  n0 = 1
  X = cbind(rep(1,dim(X)[1]),X)
  X = as.matrix(lnc.pca.df[n0:n,1:k])
  y = as.matrix(lnc.pca.df[n0:n,33])
  
  thetaGuess = runif(k*k)
  zeroGuess  = rep(0,k*k)
}
