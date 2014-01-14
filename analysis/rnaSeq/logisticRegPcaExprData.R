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
source(paste(home, "/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/exprLib.R",sep=""))
source(paste(home, "/work/research/researchProjects/encode/encode-manager/analysis/getENSGfromBiomartByRefseq.R",sep=""))
source(paste(home, "/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/pcaAnalysisWOoutlier.R",sep=""))

calcAUC <- function(prob, label){
  pre = prediction(predictions=prob, labels=label)
  per = performance(pre, "tpr", "fpr")
  AUC = (performance(pre, "auc"))@y.values[[1]] 
  AUC
}

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

#sigmoid <- function(x){ 1 / (1 + exp(-1 * x/1000))}
sigmoid <- function(x){ 1 / (1 + exp(-1 * x))}


#must accept theata as vectore
costFunction <- function(X,y,thetaVec,k){
  theta = matrix(thetaVec,k,k)
  xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
  sig = sigmoid(xTx)
  #z = -y * log(sig) - (1 - y)*log(1 - sig)
  z1 = ifelse(y == 1, -log(sig),-log(1-sig))
  z2 = ifelse(z1 == Inf, 0, z1)
  sum(z2)
}

cfLambda_old <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
    sig = abs(sigmoid(xTx))
    z1 = ifelse(y ==1, -log(sig),-log(1- sig))
    -1/length(y) * sum(ifelse(z1 == Inf, 0, z1))
  }
}


cfLambdaT <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
    sig = sigmoid(xTx)
    z1 = ifelse(y ==1, log(sig),-log(exp(xTx) + 1))
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1))
  }
}

cfLambda <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
    sig = sigmoid(xTx)
    z1 = ifelse(y ==1, log(sig),log(1- sig))
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1))
  }
}

cfLambdaReg <- function(X,y,k,lambda){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
    sig = sigmoid(xTx)
    z1 = ifelse(y ==1, log(sig),-log(exp(xTx) + 1))
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1)) + (lambda/(2 * length(y))) * sum(thetaVec ^2 )
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1)) + (lambda/(2 * length(y))) * sum(c(0,thetaVec[-1]) ^2 )
    
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
      cost = cost + tmp
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



#worksdd
gfLambdaT <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
    sig = sigmoid(xTx)
    # print(length(y),length(sig))
    s1 = sum(sig-y)
    # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
    (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k)
  } 
}

# a <- gfLambdaT(X,y,k)(thetaVec)
# gfLambdaIterate(X,y,k)(thetaVec)

gfLambdaReg <- function(X,y,k,lambda){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    regTheta = theta * (lambda / length(y))
   regTheta[1,1] = 0
    xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
    sig = sigmoid(xTx)
    # print(length(y),length(sig))
    s1 = sum(sig-y)
    # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
     (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k) + regTheta
    
  }
  
}



gfLambdaIterate <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
    sig = sigmoid(xTx)
    outMat = matrix(rep(0,k*k),k)
    for(j in 1:k){
      for(i in 1:k){ 
        tmp <- 0
        for(row in 1:length(y)){
          xxTLocal <-   X[row,i] * X[row,j] 
          tmp <- tmp + (sig[row] - y[row]) * xxTLocal 
        }
        outMat[i,j] <- tmp
      }
      
    }
    outMat * (1/length(y))
  }
}



filterXbyOptimTheta <- function(X,thetaVec,k){
  theta = matrix(thetaVec,k)
  guess = apply(X,1,function(row) t(row) %*% matrix(thetaVec,k) %*% row)
  
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
  guess = apply(X,1,function(row) t(row) %*% matrix(thetaVec,k) %*% row)
  
  # optim(fn=cfLambda(X,y,k),par=runif(k*k),method= "CG") -> o #$;o$value
  # X = as.matrix(lnc.pca.df[1:n,1:k-1])
  #  X = cbind(rep(1,dim(X)[1]),X)
  #  y = as.matrix(lnc.pca.df[1:n,33])
  
  #  guess = apply(X,1,function(row) t(row) %*% matrix(thetaVec,k) %*% row)
  #  guess
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
  n = 1500
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
  o <-optim(fn=cfLambda(X,y,k),gr=gfLambda(X,y,k),par=matrix(runif(k*k),k),method= "CG")
  # create cost function
  # create gradient
  # run optimizer
  
  results.df <- data.frame(trial=1:100,cost=rep(0,100))
  best.cost = 1000000
  best.theta =matrix(runif(k*k),k)
  for(i in 1:20){
    print(i)
    local =  optim(fn=cfLambda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS") 
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
  
  guess = apply(X.full,1,function(row) t(row) %*% matrix(best.theta,k) %*% row)
  df = data.frame(x=guess,ellipse=guess,y=y.full,label=y.full)
  # ggplot(df,aes(x=log(x),fill=factor(y)))+geom_density(alpha=I(0.6))
  df.plot= plotDistanceAwayRatio(df,label=y.full)
  df.plot$type = "ellipse"
  
  #sample(y.full,length(y.full)) 
  df.random= plotDistanceAwayRatio(df,label=sample(y.full,length(y.full)))
  df.random$type = "randomLabel"
  
  #ggplot(df.plot,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  
  df.comb = rbind(dfo,df.plot)
  ggplot(df.comb,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw() + 
    geom_smooth() + 
    facet_wrap(~type)+
    ggtitle(paste("Log Reg: rows =",n,"cols=",k,"\nPCA data used"))
  
  ggsave("~/Desktop/ellipse2.pdf",height=4,width=6)
  
  df.comb.roc = rbind(dfo,df.plot,df.random)
  ggplot(df.comb.roc,aes(x=FPR,y=TPR,color=type))+geom_line()+geom_point()+theme_bw()+geom_abline(slope=1)+
    ggtitle("ROC Curve of radius based prediction on pca\n")+xlab("False Positive Rate")+ylab("True Positive Rate")
  ggsave(file="~/Desktop/ellipse-roc1.pdf",height=4,width=5)
  
  
}

runLogReg = function(lncDf,outdir = "~/Desktop/testPCA",cols,iter=10,debug= FALSE,
                     titleMsg="",filebase=""){
  k <- length(cols) + 1
  n <- dim(lncDf)[1]
  
  if(!file.exists(outdir)){dir.create(outdir,recursive=TRUE)}
  
  
  makeOutFile <- function(x){outfile<-paste(paste(outdir,filebase,sep="/"),x,sep="");print(paste("making",outfile));gsub(" ", "",outfile)} # requires outDir & filebase
  titleWithBanner <<- function(x)paste(titleMsg,x,sep="\n")
  
  yCount <- length(which(lncDf$label == 1))
  #dataSample = c(1:94,sample(94:length(lnc.pca.df$label),(n-yCount)))
  # dataSample = 1:n
  X <- as.matrix(lncDf[,cols])
  X <- cbind(rep(1,dim(X)[1]),X)
  y <- as.matrix(lncDf[,"label"])
  
  # we need to center X:
  #  X <- apply(X,2,function(x){x - mean(x)})
  X[,1] <- 1
  
  thetaGuess <- runif(k*k)
  zeroGuess  <- rep(0,k*k)
  
  #opt.lnc = optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=thetaGuess,method= "L-BFGS-B")
  #optim(fn=cfLamda(X,y,k),gr=gfLambda(X,y,k),par=thetaGuess,method= "CG") -> o #400
  #optim(fn=cfLamda(X,y,k),par=matrix(runif(k*k),k),method= "BFGS")
  o <- optim(fn=cfLambda(X,y,k),gr=gfLambdaIterate(X,y,k),par=matrix(runif(k*k),k),method= "BFGS");o$value
  o <- optim(fn=cfLambda(X,y,k),gr=gfLambda(X,y,k),par=matrix(thetaGuess,k),method= "BFGS")
  # create cost function
  # create gradient
  # run optimizer
  
  results.df <- data.frame(trial=1:100,cost=rep(0,100))
  best.cost = 1000000
  best.theta =matrix(runif(k*k),k)
  for(i in 1:iter){
    print(i)
    local = optim(fn=cfLambda(X,y,k),par=matrix(-runif(k*k),k),method= "BFGS")
    #local =  optim(fn=cfLambda(X,y,k),par=diag(k),method= "CG") #we know this works with the getLncPcaData() fun set...see main() in this file...
    #local =  optim(fn=cfLambda(X,y,k),par=diag(k),method= "CG")
    
    results.df$trial[i] = i
    results.df$cost[i] = local$value
    if (local$value < best.cost){
      best.cost = local$value
      best.theta = local$par
    }
  }
  # guess = filterXbyOptimTheta_getGuess(X,best.theta,k)
  if("o" %in% ls() && debug == TRUE){
    best.cost = o$value
    best.theta = o$par
  }
  
  
  dfo.circle = plotDistanceAwayRatio(df=as.data.frame(X),label=y,exprCols=1:k,ellipse=FALSE)
  #ggplot(dfo,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  dfo.circle$type = "circle"
  
  guess = filterXbyOptimTheta_getGuess(X,best.theta,k)
  df.ellipse.tmp = data.frame(x=guess,ellipse=guess,y=y,label=y)
  # ggplot(df,aes(x=log(x),fill=factor(y)))+geom_density(alpha=I(0.6))
  df.ellipse = plotDistanceAwayRatio(df.ellipse.tmp,label=y,exprCols=1:k,ellipse=TRUE)
  df.ellipse$type = "ellipse"
  
  df.random.tmp  = data.frame(x=guess,ellipse=guess,y=y,label=y)
  df.random= plotDistanceAwayRatio(df.random.tmp,label=sample(y,length(y)))
  df.random$type = "randomLabel"
  
  #ggplot(df.plot,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1)+theme_bw()+geom_smooth()
  
  df.comb = rbind(dfo.circle,df.ellipse)
  ggplot(df.comb,aes(y=prec,x=sens))+geom_point()+xlim(0,1)+ylim(0,1) + 
    theme_bw() + geom_smooth() + facet_wrap(~type) +
    ggtitle(paste(titleMsg,"Precision/Sensitivity Curve for logistic regression",sep="\n"))
  ggsave(makeOutFile("elipse.pdf"),height=4,width=6)
  
  df.comb.roc = rbind(df.comb,df.random)
  ggplot(df.comb.roc,aes(x=FPR,y=TPR,color=type))+geom_line()+geom_point()+theme_bw()+geom_abline(slope=1)+
    ggtitle(paste(titleMsg,"ROC Curve of radius based prediction on pca\n", sep="\n"))+xlab("False Positive Rate")+ylab("True Positive Rate")
  ggsave(makeOutFile("elipse-ROC.pdf"),height=4,width=5)
  
  
}

runLogRegTheta = function(X,k,y,cols,iter=10,reg= FALSE,lambda = 1){
  if (TRUE == reg){
    costFunction <- cfLambdaReg(X,y,k,lambda)
    
  }
  else {  
    costFunction <- cfLambda(X,y,k)
  }
  
  thetaGuess <- runif(k*k)
  zeroGuess  <- rep(0,k*k)
  
  results.df <- data.frame(trial=1:100,cost=rep(0,100))
  best.cost = 1000000
  best.theta =matrix(runif(k*k),k)
  for(i in 1:iter){
    print(i)
    local = optim(fn=costFunction,par=matrix(-runif(k*k),k),method= "BFGS")
    results.df$trial[i] = i
    results.df$cost[i] = local$value
    if (local$value < best.cost){
      best.cost = local$value
      best.theta = local$par
    }
  }
  # guess = filterXbyOptimTheta_getGuess(X,best.theta,k)
  best.theta;
}

testTheta <- function(X,k,y,cols, theta){
  
  predictLabel <- ifelse(apply(X,1,function(x){t(x) %*% theta %*% x} ) < 0, 0, 1)
  predictY <- ifelse(predictLabel == y, 1, 0)
  sum(predictY)/ length(predictY)
}

testThetaTPR <- function(X,k,y,cols, theta){
  
  predictLabel <- ifelse(apply(X,1,function(x){t(x) %*% theta %*% x} ) < 0, 0, 1)
  predictY <- ifelse(predictLabel == y, 1, 0)
  sum(predictY[which(y == 1)])/ length(predictY[which(y == 1)])
}

testThetaStats <- function(X,k,y,cols, theta){
  probs <- sigmoid(apply(X,1,function(x){t(x) %*% theta %*% x} ))
  AUC <- calcAUC(probs, y)
  predictLabel <- ifelse(apply(X,1,function(x){t(x) %*% theta %*% x} ) < 0, 0, 1)
  predictY <- ifelse(predictLabel == y, 1, 0)
  TP <- length(which(predictY == 1 & y == 1))
  TN <- length(which(predictY == 1 & y == 0))
  FN <- length(which(predictY == 0 & y == 1))
  FP <- length(which(predictY == 0 & y == 0))
  list(TP=TP,TN=TN,FN=FN,FP=FP,AUC=AUC)
  
}

trainAndTestLogReg <- function(lncDf, ratio, cols, mainEffects=TRUE,reg=FALSE,lambda=1,iter=1,glmTypeOutput=FALSE){
  
  totalCount = dim(lncDf)[1]
  
  k <- length(cols) + 1
  n <- dim(lncDf)[1]
  yCount <- length(which(lncDf$label == 1))
  
  X <- as.matrix(lncDf[,cols])
  X <- cbind(rep(1,dim(X)[1]),X)
  y <- as.matrix(lncDf[,"label"])
  
  # we need to center X:
  #  X <- apply(X,2,function(x){x - mean(x)})
  
  if (TRUE == mainEffects){
    X[,1] <- 1
  } else{
    X <- X[,-1]
    k <- k - 1
  }
  
  #scale X
  X <- X / max(X)
  
  train = sample(1:totalCount, totalCount * ratio)
  #theta = runLogRegTheta(X[train,],k,y,cols,iter=iter,reg=reg,lambda=lambda)
  theta = matrix(runNLM(X[train,],y[train],k,reg,lambda)$estimate,k)
  
  if(glmTypeOutput){
    probs <- sigmoid(apply(X[-train,],1,function(x){t(x) %*% theta %*% x} ))
    getStatsFromGlmModel(probs, y[-train], knn=FALSE)
    
  } else {
  

  
  #testThetaStatsEllipse(X,k,y,cols, theta)
  
  four <- testThetaStats(X,k,y,cols, theta) # TP=four$TP,TN=four$TN,FP=four$FP,FN=four$FN
  list(R=testTheta(X[-train,],k,y[-train],cols,theta),TPR =testThetaTPR(X[-train,],k,y[-train],cols,theta),TP=four$TP,TN=four$TN,FP=four$FP,FN=four$FN)
  }
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


test1 = function(outdir = "/home/wespisea/work/research/researchProjects/encode/encode-manager/plots/fullAnalysisExperiment/test/logReg/sampleData/"){
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  
  nn = 1000
  lower.df= data.frame(x0 = rep(1,nn),
                       x1 = runif(nn)*2 -1 ,
                       x2 = runif(nn)*2 - 1,
                       y = rep(0,nn))
  lower.df$dist <- apply(lower.df, 1, function(x)sqrt(x[["x2"]]^2  + x[["x1"]]^2))
  lower.df <- lower.df[which(lower.df$dist < 1),]
  lowerMax =  ifelse(dim(lower.df)[1] > 100, 100, dim(lower.df)[1])
  upper.df= data.frame(x0 = rep(1,nn),
                       x1 = runif(nn)*3-1.5,
                       x2 = runif(nn)*3-1.5,
                       y = rep(1,nn))
  upper.df$dist <- apply(upper.df, 1, function(x)sqrt(x[["x2"]]^2  + x[["x1"]]^2))
  upper.df <- upper.df[which(upper.df$dist < 1.5 & upper.df$dist > 1),]
  upperMax =  ifelse(dim(upper.df)[1] > 100, 100, dim(upper.df)[1])
  
  X.df = rbind(lower.df[1:lowerMax,],upper.df[1:upperMax,])
  X.test = as.matrix(X.df[,1:3])
  y.test = X.df$y
  k.test = 3
  
  train = sample(seq_along(X.df$y),0.7 *length(X.df$y))
  
  ggplot(X.df,aes(x1,x2,color=factor(y))) + geom_point() + theme_bw()+ggtitle("data")
  ggsave(paste(outdir,"dataPlot.pdf"))
  
  ggplot(X.df[train,],aes(x1,x2,color=factor(y))) + geom_point() + 
    theme_bw()+ggtitle("Training data")
  ggsave(paste(outdir,"dataPlot-train.pdf"))
  
  ggplot(X.df[-train,],aes(x1,x2,color=factor(y))) + geom_point() + 
    theme_bw()+ggtitle("Testing data")
  ggsave(paste(outdir,"dataPlot-test.pdf"))
  
  theta = matrix(runNLM(X[train,],y.test[train],k=3,reg=TRUE,0)$estimate,3)
  X.df$predict <- ifelse(apply(X,1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  X.df$correct <- ifelse(X.df$predict == X.df$y, 1, 0)       
  #t1 = data.frame(x=(1:1000-500)/10)
  #t1$y = apply(t1,1,function(c)test1(diag(k.test),X.test,k.test,c))
  
  trainR = sum(X.df[train,"correct"])/dim(X.df[train,])[1]
  ggplot(X.df[train,],aes(x1,x2,color=factor(y),size=1-correct)) +
    geom_point() + theme_bw()+ ggtitle(paste("test performance -- wrong points are large\nR=",trainR,sep="")) 
  ggsave(paste(outdir,"dataPlot-test-Incorrect.pdf"))
  
  testR = sum(X.df[-train,"correct"])/dim(X.df[-train,])[1]
  ggplot(X.df[-train,],aes(x1,x2,color=factor(y),size=1-correct)) +
    geom_point() + theme_bw()+ ggtitle(paste("test performance -- wrong points are large\nR=",testR,sep="")) 
  ggsave(paste(outdir,"dataPlot-test-Incorrect.pdf"))
  
  analysis.df <- data.frame(expand.grid(x1=seq(-2,2,by=0.1),x2=seq(-2,2,by=0.1)))
  analysis.df$x0 <- 1
  theta <- matrix(runNLM(X[train,],y.test[train],k=3,reg=TRUE,0)$estimate,3)
  analysis.df$predict <- ifelse(apply(as.matrix(analysis.df[c("x0","x1","x2")],),1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  ggplot(analysis.df, aes(x1,x2,fill=factor(predict)))+geom_tile() + theme_bw()+
    ggtitle("Decision boundary for training\nLambda = 0")
  ggsave(paste(outdir,"dataPlot-decisionBoundary-lambda=0.pdf"))
  
  analysis.df <- data.frame(expand.grid(x1=seq(-2,2,by=0.1),x2=seq(-2,2,by=0.1)))
  analysis.df$x0 <- 1
  theta <- matrix(runNLM(X[train,],y.test[train],k=3,reg=TRUE,10)$estimate,3)
  analysis.df$predict <- ifelse(apply(as.matrix(analysis.df[c("x0","x1","x2")],),1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  ggplot(analysis.df, aes(x1,x2,fill=factor(predict)))+geom_tile() + theme_bw()+
    ggtitle("Decision boundary for training\nLambda = 10")
  ggsave(paste(outdir,"dataPlot-decisionBoundary-lambda=10.pdf"))
  
  analysis.df <- data.frame(expand.grid(x1=seq(-2,2,by=0.1),x2=seq(-2,2,by=0.1)))
  analysis.df$x0 <- 1
  theta <- matrix(runNLM(X[train,],y.test[train],k=3,reg=TRUE,50)$estimate,3)
  analysis.df$predict <- ifelse(apply(as.matrix(analysis.df[c("x0","x1","x2")],),1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  ggplot(analysis.df, aes(x1,x2,fill=factor(predict)))+geom_tile() + theme_bw()+
    ggtitle("Decision boundary for training\nLambda = 50")
  ggsave(paste(outdir,"dataPlot-decisionBoundary-lambda=50.pdf"))
  
}

test2 = function(outdir = "/home/wespisea/work/research/researchProjects/encode/encode-manager/plots/fullAnalysisExperiment/test/logReg/sampleData2/"){
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  xval = 1
  nn = 1000
  lower.df= data.frame(x0 = rep(xval,nn),
                       x1 = runif(nn)*0.5 -0.25 ,
                       x2 = runif(nn)*0.5 - 0.25,
                       y = rep(0,nn))
  lower.df$dist <- apply(lower.df, 1, function(x)sqrt(x[["x2"]]^2  + x[["x1"]]^2))
  lower.df <- lower.df[which(lower.df$dist < 1),]
  lowerMax =  ifelse(dim(lower.df)[1] > 200, 200, dim(lower.df)[1])
  ones = 200
  upper.df= data.frame(x0 = rep(xval,ones),
                       x1 = runif(ones)*0.5,
                       x2 = runif(ones)*0.5,
                       y = rep(1,ones))
  upper.df[1:50, "x2"] <- upper.df[1:50, "x2"] + 1
  upper.df[1:50, "x1"] <- upper.df[1:50, "x1"] - 0.25
  upper.df[51:100, "x2"] <- upper.df[51:100, "x2"] - 1.5
  upper.df[51:100, "x1"] <- upper.df[51:100, "x1"] - 0.25
  upper.df[101:150, "x1"] <- upper.df[101:150, "x1"] + 1
  upper.df[101:150, "x2"] <- upper.df[101:150, "x2"] - 0.25
  upper.df[151:200, "x1"] <- upper.df[151:200, "x1"] - 1.5
  upper.df[151:200, "x2"] <- upper.df[151:200, "x2"] - 0.25
  upper.df$dist = 0 
  
  X.df = rbind(lower.df[1:200,], upper.df)
  X.test = as.matrix(X.df[,c("x0","x1","x2")])
  y.test = X.df$y
  k.test = 3
  
  train = 1:300
  
  ggplot(X.df,aes(x1,x2,color=factor(y))) + geom_point() + theme_bw()+ggtitle("data")+ xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot.pdf"))
  
  ggplot(X.df[train,],aes(x1,x2,color=factor(y))) + geom_point() + 
    theme_bw()+ggtitle("Training data") + xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-train.pdf"))
  
  ggplot(X.df[-train,],aes(x1,x2,color=factor(y))) + geom_point() + 
    theme_bw()+ggtitle("Testing data")+ xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-test.pdf"))
  
  theta = matrix(runNLM(X.test[train,],y.test[train],k=3,reg=TRUE,0)$estimate,3)
  X.df$predict <- ifelse(apply(X.test,1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  X.df$correct <- ifelse(X.df$predict == X.df$y, 1, 0)       
  #t1 = data.frame(x=(1:1000-500)/10)
  #t1$y = apply(t1,1,function(c)test1(diag(k.test),X.test,k.test,c))
  
  trainR = sum(X.df[train,"correct"])/dim(X.df[train,])[1]
  ggplot(X.df[train,],aes(x1,x2,color=factor(y),size=1-correct)) +
    geom_point() + theme_bw()+ ggtitle(paste("test performance -- wrong points are large\nR=",trainR,sep="")) + xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-train-Incorrect.pdf"))
  
  testR = sum(X.df[-train,"correct"])/dim(X.df[-train,])[1]
  ggplot(X.df[-train,],aes(x1,x2,color=factor(y),size=1-correct)) +
    geom_point() + theme_bw()+ ggtitle(paste("test performance -- wrong points are large\nR=",testR,sep="")) + xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-test-Incorrect.pdf"))
  
  analysis.df <- data.frame(expand.grid(x1=seq(-2,2,by=0.1),x2=seq(-2,2,by=0.05)))
  analysis.df$x0 <- xval
  theta <- matrix(runNLM(X.test[train,],y.test[train],k=3,reg=TRUE,0)$estimate,3)
  analysis.df$predict <- ifelse(apply(as.matrix(analysis.df[c("x0","x1","x2")],),1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  ggplot(analysis.df, aes(x1,x2,fill=factor(predict)))+geom_tile() + theme_bw()+
    ggtitle("Decision boundary for training\nLambda = 0")+ xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-decisionBoundary-lambda=0.pdf"))
  
  analysis.df <- data.frame(expand.grid(x1=seq(-2,2,by=0.1),x2=seq(-2,2,by=0.05)))
  analysis.df$x0 <- xval
  analysis.df$dist <-  apply(analysis.df, 1, function(x)sqrt(x[["x2"]]^2  + x[["x1"]]^2))
  theta <- matrix(runNLM(X.test[train,],y.test[train],k=3,reg=TRUE,10)$estimate,3)
  analysis.df$predict <- ifelse(apply(as.matrix(analysis.df[c("x0","x1","x2")],),1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  analysis.df$val <- apply(as.matrix(analysis.df[c("x0","x1","x2")],),1,function(x)t(x) %*% matrix(theta,3) %*% x)
  ggplot(analysis.df, aes(x1,x2,fill=factor(predict)))+geom_tile() + theme_bw()+
    ggtitle("Decision boundary for training\nLambda = 10")+ xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-decisionBoundary-lambda=10.pdf"))
  
  analysis.df <- data.frame(expand.grid(x1=seq(-10,10,by=0.1),x2=seq(-10,10,by=0.05)))
  analysis.df$x0 <- xval
  theta <- matrix(runNLM(X.test[train,],y.test[train],k=3,reg=TRUE,50)$estimate,3)
  analysis.df$predict <- ifelse(apply(analysis.df[c("x0","x1","x2")],1,function(x)t(x) %*% matrix(theta,3) %*% x) < 0, 0,1)
  ggplot(analysis.df, aes(x1,x2,fill=factor(predict)))+geom_tile() + theme_bw()+
    ggtitle("Decision boundary for training\nLambda = 50")+ xlim(-2,2)+ ylim(-2,2)
  ggsave(paste(outdir,"dataPlot-decisionBoundary-lambda=50.pdf"))
  
}


runNLM <- function(X,y,k,reg,lambda,useGradient=TRUE){  
  if (TRUE == reg){
    dim(X)
    a <- cfLambdaReg(X,y,k,lambda)
    b <- gfLambdaReg(X,y,k,lambda)
  } else {
    a <- cfLambdaT(X,y,k)
    b <- gfLambdaT(X,y,k)
  }
  if(useGradient){
  fgh <- function(x){
    res <- a(x)
    attr(res, "gradient") <- b(x)
    return(res)
  }} else {
    fgh <- function(x){
      res <- a(x)
      #attr(res, "gradient") <- b(x)
      return(res)
    }
    
    
  }
  nlm(f= fgh, p = runif(k*k),check.analyticals=TRUE)
  
}

nlmTest <- function(){
  # X <- X / max(X)
  a <- cfLambdaT(X,y,k)
  yy <- ifelse(y == 1, 0,1)
  b <- gfLambdaT(X,y,k)
  fgh <- function(x){
    res <- a(x)
    attr(res, "gradient") <- b(x)
    return(res)
  }
  #thetaVec = runif(k*k)
  nlm(f= fgh, p = runif(k*k),check.analyticals=TRUE)
  
}
nlmTestReg <- function(){
  # X <- X / max(X)
  a <- cfLambdaReg(X,y,k,lambda=1000)
  b <- gfLambdaReg(X, ifelse(y == 1, 1,1),k,lambda=1000)
  fgh <- function(x){
    res <- a(x)
    attr(res, "gradient") <- b(x)
    return(res)
  }
  #thetaVec = runif(k*k)
  nlm(f= fgh, p = runif(k*k),check.analyticals=TRUE)
}
