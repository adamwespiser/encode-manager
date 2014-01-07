home <- Sys.getenv("HOME")

projectDir <- paste(home,"/work/research/researchProjects/encode/encode-manager",sep="")
setwd(projectDir)
getwd()

#boiler plate helpers
getFullPath <- function(projectPath){ paste(projectDir,projectPath,sep="/") }
exportAsTable <- function(df, file){ write.table(df,file=file,quote=FALSE, row.names=FALSE,sep="\t") }
clear <- function(save.vec=c()){ ls.vec <- ls(globalenv());del.vec <-setdiff(ls.vec,c(save.vec,"clear")); rm(list=del.vec,pos=globalenv())}
readInTable <- function(file) read.table(file=file,stringsAsFactors=FALSE,header=TRUE)

sigmoid <- function(x){ 1 / (1 + exp(-1 * x))}
library(class)
library(MASS)
library(LiblineaR)
library(ROCR)

cfLambdaT <- function(X,y,k){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
    sig = sigmoid(xTx)
    z1 = ifelse(y ==1, log(sig),-log(exp(xTx) + 1))
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
  }
}
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
gfLambdaReg <- function(X,y,k,lambda){
  function(thetaVec){
    theta = matrix(thetaVec,k)
    regTheta = theta * (lambda / length(y))
    xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
    sig = sigmoid(xTx)
    # print(length(y),length(sig))
    s1 = sum(sig-y)
    # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
    (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k) + regTheta
  }
  
}

cfLambdaTStraightVec <- function(X,y){
  function(thetaVec){
    
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    sig = sigmoid(xT)
    z1 = ifelse(y ==1, log(sig),-log(exp(xT) + 1))
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1))
  }
}
cfLambdaRegStraightVec <- function(X,y,lambda){
  function(thetaVec){
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    sig = sigmoid(xT)
    z1 = ifelse(y ==1, log(sig),-log(exp(xT) + 1))
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1)) + (lambda/(2 * length(y))) * sum(thetaVec ^2 )
  }
}
gfLambdaTStraightVec <- function(X,y){
  function(thetaVec){
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    sig = sigmoid(xT)
    # print(length(y),length(sig))
    s1 = sum(sig-y)
    # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
   # (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k)
    (1/length(y)) * s1 * thetaVec
  } 
}
gfLambdaRegStraightVec <- function(X,y,lambda){
  function(thetaVec){
    regTheta = thetaVec * (lambda / length(y))
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    sig = sigmoid(xT)
    # print(length(y),length(sig))
    s1 = sum(sig-y)
    # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
    #(1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k) + regTheta
    (1/length(y)) * s1 * thetaVec + regTheta
  }
  
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

runNlmTest <- function(X,y,grad=TRUE){
  # X <- X / max(X)
  k <- dim(X)[2]
  a <- cfLambdaTStraightVec(X,y)
  b <- gfLambdaTStraightVec(X,y)
  if(grad){
    fgh <- function(x){
      res <- a(x)
      attr(res, "gradient") <- b(x)
      return(res)
    }
  } else {
    fgh <- function(x){
      res <- a(x)
      # attr(res, "gradient") <- b(x)
      return(res)
    }
  }
  #thetaVec = runif(k*k)
  nlm(f= fgh, p = runif(k),check.analyticals=TRUE)
  
}

runNlmTestReg <- function(X,y,lambda,grad=TRUE){
  # X <- X / max(X)
  k <- dim(X)[2]
  a <- cfLambdaRegStraightVec(X,y,lambda)
  b <- gfLambdaRegStraightVec(X,y,lambda)
  if(grad){
    fgh <- function(x){
      res <- a(x)
      attr(res, "gradient") <- b(x)
      return(res)
    }
  } else {
    fgh <- function(x){
      res <- a(x)
      # attr(res, "gradient") <- b(x)
      return(res)
    }
  }
  #thetaVec = runif(k*k)
  nlm(f= fgh, p = runif(k),check.analyticals=TRUE)
  
}


mapFeature <- function(X,degree=6){
  if (dim(X)[2] != 2){
    warning("supplied matrix has more then 2 columns...\nreturning X")
    return(X)
  }
  for (i in 1:degree){
    for (j in 0:i){
      X <- cbind(X, (X[,1]^ (i-j)) * (X[,2]^ j) ) 
      
    }    
  }
  X[,-(1:2)]
}

testThetaStats <- function(X, y, thetaVec){
  k <- dim(X)[2]
  probs <- sigmoid(apply(X,1,function(row)t(row)%*% thetaVec))
  AUC <- calcAUC(probs, y)
  predictLabel <- ifelse(apply(X,1,function(x){t(x) %*% thetaVec} ) < 0, 0, 1)
  predictY <- ifelse(predictLabel == y, 1, 0)
  TP <- length(which(predictY == 1 & y == 1))
  TN <- length(which(predictY == 1 & y == 0))
  FN <- length(which(predictY == 0 & y == 1))
  FP <- length(which(predictY == 0 & y == 0))
  sens = TP / (TP + FN)
  prec = TP / (TP + FP)
  F1 = 2*(prec*sens)/(prec+sens)
  list(TP=TP,TN=TN,FN=FN,FP=FP,AUC=AUC,sens=sens,prec=prec,F1=F1)
  
}



predictWithThetaOnX <- function(X,thetaVec){
  probs <- apply(X,1,function(row)t(row)%*% thetaVec)
  probs
}


testWithMlClassData <- function(){
  mldata <- read.csv(getFullPath("data/logRegTest/mlclass-ex2/ex2data2.txt"),sep=",",header=FALSE)
  X <- cbind(rep(1,dim(mldata)[1]), mldata[,1:2])
  y <- mldata[,3]

  nlm.out <- runNlmTest(X,ifelse(y==1,0,1))
  predictWithThetaOnX(X,nlm.out$estimate)
}


plotPredictionML <- function(inputFile=getFullPath("data/logRegTest/mlclass-ex2/ex2data1.txt"),
                                     outputFile=getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data1.pdf"),
                                     stepSize=0.02,
                                     degree=1,exprName="ex2data1",
                                     useReg=FALSE,
                                     lambda=1,
                                     useGradient=TRUE,ySwap=TRUE,scaleData=FALSE){
  mldata <- read.csv(inputFile,sep=",",header=FALSE)
  
  if(scaleData){
  X.scale <- data.frame(scale(mldata[,1:2],center=TRUE,scale=TRUE))
  } else {
  X.scale <- mldata[,1:2]
  }
  
  if(ySwap){
    y <- ifelse(y==1,0,1)
  } else {
    y
  }
  
  X.feature <- mapFeature(X.scale,degree)
  X <- cbind(rep(1,dim(mldata)[1]), X.feature)
  #X <- mldata[,1:2]
   
  y <- mldata[,3]
  
  if(!useReg){
    nlm.out <- runNlmTest(X,y,grad=useGradient)
  } else {
    print(paste("use reg lambda=",lambda,"\n",sep=""))
    nlm.out <- runNlmTestReg(X,y,grad=useGradient,lambda=lambda)
  }  
  thetaVec <- nlm.out$estimate
  
  df <- data.frame(cbind(X,y))
  col.names <-paste0("X",seq_along(names(X))-1)
  colnames(df) <- c(col.names,"y")
  df$logLik <- apply(df[col.names],1,function(row)t(row) %*% thetaVec)
  df$predict <- ifelse(df$logLik < 0, 0, 1 )
  stats.list <- testThetaStats(df[col.names], y, thetaVec)
  # create grid for the prediction
  grid.df <- data.frame(expand.grid(X1=seq(min(df$X1),max(df$X1),by=stepSize),X2=seq(min(df$X2),max(df$X2),by=stepSize)))
  grid.mat <- mapFeature(as.matrix(grid.df),degree=degree)
  grid.df <- as.data.frame(cbind(rep(1,dim(grid.df)[1]), grid.mat))
  colnames(grid.df) <- col.names
  grid.df$logLik <- apply(grid.df[col.names],1,function(row)t(row) %*% thetaVec)
  grid.df$predict <- ifelse(grid.df$logLik < 0, 0, 1 )
  
  # figure out the title, with stats, dataset info, and lambda info(if regularized)
  stats.title <- with(stats.list,paste("\nAUC=",round(AUC,digits=2)," F1=",round(F1,digits=2),
                                       " sens=",round(sens,digits=2)," prec=",round(prec,digits=2),"\n",sep=""))
  title <- paste("mlclass ex 2 training dataset=",exprName," poly degree=",degree,stats.title,"\ncost fn = ",nlm.out$minimum,sep="")
  title <- ifelse(useReg, paste(title, "\nreg'd w/ lambda =",lambda,sep=""),title)
  
  ggplot(grid.df, aes(x= X1,y= X2,fill=factor( predict)))+geom_raster(alpha=I(0.4))+
    theme_bw()+
    layer( data = df,
           mapping = aes(x= X1,y= X2,color=factor( y),size=3),
           geom="point") +
    ggtitle(title)
  
  ggsave(file=outputFile,height=6,width=6)
}


runTests <- function(){
#   plotPredictionFirstOrder(inputFile=getFullPath("data/logRegTest/mlclass-ex2/ex2data1.txt"),
#                           outputFile=getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data1_degree=1.pdf"),
#                            degree=1,exprName=)
#   
#   plotPredictionFirstOrder(inputFile=getFullPath("data/logRegTest/mlclass-ex2/ex2data2.txt"),
#                            outputFile=getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data2_degree=1.pdf"),
#                            degree=1)
#   
  for (dataset in 2:2){
    dataset.file <- paste(getFullPath("data/logRegTest/mlclass-ex2/ex2data"), dataset, ".txt",sep="")
    for (i in 1:6){
      output.file <- paste(getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data"), dataset,"_degree=",i, ".pdf",sep="")
      plotPredictionML(inputFile=dataset.file,
                               outputFile=output.file,
                               degree=i,
                               exprName = paste("ex2data",i,sep=""),useGradient=TRUE)
    }
    
  }
  
  lambda.vals <- c(0,1,100)
  for(l in lambda.vals){
    dataset.file <- paste(getFullPath("data/logRegTest/mlclass-ex2/ex2data"), 2, ".txt",sep="")
    output.file <- paste(getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data2"),"_lambda=",l,"_degree=6",".pdf",sep="")
    plotPredictionML(inputFile=dataset.file,
                     outputFile=output.file,
                     degree=6,
                     exprName = paste("ex2data",i,sep=""),useGradient=FALSE,
                     useReg = TRUE, lambda = l)
  }
  
}






