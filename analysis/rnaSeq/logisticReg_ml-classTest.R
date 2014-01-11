home <- Sys.getenv("HOME")
source(paste(home,"/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/logisticRegPcaExprData.R",sep=""))

########### function to use is: plotPredictionML() 

############### to run analysis in github:
## rm(list=unlist(ls()));source('~/work/research/researchProjects/encode/encode-manager/analysis/rnaSeq/logisticReg_ml-classTest.R');runTests()

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

# cfLambdaT <- function(X,y,k){
#   function(thetaVec){
#     theta = matrix(thetaVec,k)
#     xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
#     sig = sigmoid(xTx)
#     z1 = ifelse(y ==1, log(sig),-log(exp(xTx) + 1))
#     -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1))
#   }
# }
# cfLambdaReg <- function(X,y,k,lambda){
#   function(thetaVec){
#     theta = matrix(thetaVec,k)
#     xTx = apply(X,1,function(row){t(row) %*% theta %*% row})
#     sig = sigmoid(xTx)
#     z1 = ifelse(y ==1, log(sig),-log(exp(xTx) + 1))
#     -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1)) + (lambda/(2 * length(y))) * sum(thetaVec ^2 )
#   }
# }
# gfLambdaT <- function(X,y,k){
#   function(thetaVec){
#     theta = matrix(thetaVec,k)
#     xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
#     sig = sigmoid(xTx)
#     # print(length(y),length(sig))
#     s1 = sum(sig-y)
#     # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
#     (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k)
#   } 
# }
# gfLambdaReg <- function(X,y,k,lambda){
#   function(thetaVec){
#     theta = matrix(thetaVec,k)
#     regTheta = theta * (lambda / length(y))
#     xTx = apply(X,1,function(row)t(row) %*% theta %*% row)
#     sig = sigmoid(xTx)
#     # print(length(y),length(sig))
#     s1 = sum(sig-y)
#     # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
#     (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k) + regTheta
#   }
#   
# }

cfLambdaTStraightVec <- function(X,y){
  function(thetaVec){
    X <- as.matrix(X)
    y <- as.matrix(y)
    
    thetaVec <- as.matrix(thetaVec)
    h <- sigmoid(X %*% thetaVec)
    m <- length(y)
    
    
    (t(-y)%*%log(h) - t(1-y)%*%log(1-h))/m ;
  }
}
cfLambdaTStraightVec_old <- function(X,y){
  function(thetaVec){
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    sig = sigmoid(xT)
    z1 = ifelse(y == 1, log(sig),-log(exp(xT) + 1))
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1))
  }
}


cfLambdaRegStraightVec <- function(X,y,lambda){
  function(thetaVec){
   # thetaVec[1] <- 0
    X <- as.matrix(X)
    y <- as.matrix(y)
    
    thetaVec <- as.matrix(thetaVec)
    h <- sigmoid(X %*% thetaVec)
    m <- length(y)
    
    #p = lambda*(t(thetaVec)%*%thetaVec)/(2*m);
    p = lambda*(sum(thetaVec[-1] ^ 2))/(2*m)
  (t(-y)%*%log(h) - t(1-y)%*%log(1-h))/m + p;


  }
}


cfLambdaRegStraightVec_old <- function(X,y,lambda){
  function(thetaVec){
    thetaVec[1] <- 0
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    h = sigmoid(xT)
    z1 = ifelse(y == 1, log(h),-log(exp(xT) + 1))
    reg <- (lambda * sum(thetaVec ^2 ) /(2 * length(y))) 
    -1/length(y) * sum(ifelse(z1 == -Inf, 0, z1)) + reg
  }
}

gfLambdaTStraightVec <- function(X,y){
  function(thetaVec){
    #thetaVec[1] <- 0
    X <- as.matrix(X)
    y <- as.matrix(y)
    
    thetaVec <- as.matrix(thetaVec)
    h <- sigmoid(X %*% thetaVec)
    (t(X)%*%(h - y))/length(y);
  } 
}


gfLambdaTStraightVec_old <- function(X,y){
  function(thetaVec){
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    h = sigmoid(xT)
    # print(length(y),length(sig))
    s1 = sum(h-y)
    # as.vector(unlist(numcolwise(mean)(as.data.frame(t(apply(X,1,function(x){x %*% t(x)})))) * dim(X)[1])) * s1
   # (1/length(y)) * matrix(apply(X,1,function(x){as.vector(tcrossprod(x))}) %*% (sig-y),k)
    (1/length(y)) * s1 * thetaVec
  } 
}
gfLambdaRegStraightVec <- function(X,y,lambda){
  function(thetaVec){
 # thetaVec[1] <- 0
  X <- as.matrix(X)
  y <- as.matrix(y)
  
  thetaVec <- as.matrix(thetaVec)
  h <- sigmoid(X %*% thetaVec)
  (t(X)%*%(h - y)+c(0,lambda*thetaVec[-1]))/length(y);
  } 
}
gfLambdaRegStraightVec_old <- function(X,y,lambda){
  function(thetaVec){
    regTheta = thetaVec * (lambda / length(y))
    xT = apply(X,1,function(row)t(row)%*% thetaVec)
    sig = sigmoid(xT)
    #print(length(y),length(sig))
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

runNlmTestML <- function(X,y,grad=TRUE){
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
  nlm(f= fgh, p = runif(k)*0,check.analyticals=TRUE)
  
}

runNlmTestRegML <- function(X,y,lambda,grad=TRUE){
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
  nlm(f= fgh, p = runif(k)*0,check.analyticals=TRUE)
  
}

runNLMChooseReg <- function(X,y,grad=TRUE,lambda=1,useReg=TRUE){
  if (useReg){
    runNlmTestRegML(X,y,lambda,grad=grad)
  } else {
    runNlmTestML(X,y,grad=grad)
  }
  
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
testThetaStatsEllipse <- function(X, y, thetaVec){
  k <- dim(X)[2]
  probs <- sigmoid(apply(X,1,function(row)t(row)%*% matrix(thetaVec,3) %*% row))
  AUC <- calcAUC(probs, y)
  predictLabel <- ifelse(apply(X,1,function(row)t(row)%*% matrix(thetaVec,3) %*% row) < 0, 0, 1)
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

plotPredictionMLEllipseTransform <- function(inputFile=getFullPath("data/logRegTest/mlclass-ex2/ex2data1.txt"),
                             outputFile=getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data1.pdf"),
                             stepSize=0.02,
                             degree=1,exprName="ex2data1",
                             useReg=FALSE,
                             lambda=1,
                             useGradient=TRUE,ySwap=TRUE,scaleData=FALSE,useEllipse=FALSE,msg){
  mldata <- read.csv(inputFile,sep=",",header=FALSE)
  
  if(scaleData){
    X.scale <- data.frame(scale(mldata[,1:2],center=TRUE,scale=TRUE))
  } else {
    X.scale <- mldata[,1:2]
  }
  
  if(ySwap){
    y <- ifelse(y==1,0,1)
  } 
  
  #X.feature <- mapFeature(X.scale,degree)
  X <- cbind(rep(1,dim(mldata)[1]), X.scale)
  X <- colsTransformEllipse(X,k=3)
  dim(X)
  #X <- mldata[,1:2]
  
  y <- mldata[,3]
  if(!useReg){
    nlm.out <- runNlmTestML(X,y,grad=useGradient)
  } else {
    print(paste("use reg lambda=",lambda,"\n",sep=""))
    nlm.out <- runNlmTestRegML(X,y,grad=useGradient,lambda=lambda)
  }  
  df <- data.frame(cbind(X,y))
  thetaVec <- nlm.out$estimate
  
  
  col.names <-paste0("X",seq_along(names(X))-1)
  colnames(df) <- c(col.names,"y")

  df$logLik <- apply(X,1,function(row)t(row) %*% thetaVec)
    
  stats.list <- testThetaStats(df[col.names], y, thetaVec)
  df$predict <- ifelse(df$logLik < 0, 0, 1 )
  # create grid for the prediction
  grid.df <- data.frame(expand.grid(X1=seq(min(df$X1),max(df$X1),by=stepSize),X2=seq(min(df$X2),max(df$X2),by=stepSize)))
  grid.df <- cbind(rep(1,dim(grid.df)[1]),grid.df)
  grid.df <- colsTransformEllipse(grid.df,k=3)
  
  colnames(grid.df) <- col.names
  grid.df$logLik <- apply(grid.df,1,function(row)t(row) %*% thetaVec)
  grid.df$predict <- ifelse(grid.df$logLik < 0, 0, 1 )
  
  # figure out the title, with stats, dataset info, and lambda info(if regularized)
  stats.title <- with(stats.list,paste("\nAUC=",round(AUC,digits=2)," F1=",round(F1,digits=2),
                                       " sens=",round(sens,digits=2)," prec=",round(prec,digits=2),"\n",sep=""))
  title <- paste("mlclass ex 2 training dataset=",exprName," poly degree=",degree,stats.title,"cost fn = ",round(nlm.out$minimum,digits=2),sep="")
  title <- ifelse(useReg, paste(title, "\nreg'd w/ lambda =",lambda,sep=""),title)
  title <- paste(title,msg,sep="")
  ggplot(grid.df, aes(x= X1,y= X2,fill=factor( predict)))+geom_raster(alpha=I(0.4))+
    theme_bw()+
    layer( data = df,
           mapping = aes(x= X1,y= X2,color=factor( y),size=3),
           geom="point") +
    ggtitle(title)
  
  ggsave(file=outputFile,height=8,width=8)
}
plotPredictionML <- function(inputFile=getFullPath("data/logRegTest/mlclass-ex2/ex2data1.txt"),
                             outputFile=getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data1.pdf"),
                             stepSize=0.02,
                             degree=1,exprName="ex2data1",
                             useReg=FALSE,
                             lambda=1,
                             useGradient=TRUE,ySwap=FALSE,scaleData=FALSE,
                             useEllipse=FALSE,ellipseTransform=FALSE,msg){
  mldata <- read.csv(inputFile,sep=",",header=FALSE)
  if(useEllipse && ellipseTransform){
    warning("useEllipse and ellipseTransform set to true -> pick one or the other")
    return(-1)
  }
  X.scale <- mldata[,1:2]
  y <- mldata[,3]
  
  if(scaleData){
    X.scale <- data.frame(X.scale,center=TRUE,scale=TRUE)
  } 
  
  if(ySwap){y <- ifelse(y==1,0,1)} 
  
  if(ellipseTransform){
    X <- cbind(rep(1,dim(mldata)[1]), X.scale)
    X <- colsTransformEllipse(X,k=3)
    
  } else {
    X.feature <- mapFeature(X.scale,degree)
    X <- cbind(rep(1,dim(mldata)[1]), X.feature)
  }
  
  
  
  df <- data.frame(cbind(X,y))

  if(useEllipse && !ellipseTransform){
    X.ellipse <- X.scale[,1:2]
    X.ellipse <- cbind(rep(1,dim(X.ellipse)[1]),X.ellipse)
    df <- data.frame(cbind(X,y))
    nlm.out <- runNLM(X=as.matrix(X.ellipse),y=y,k=3,reg=useReg,useGradient=useGradient,lambda=lambda)
  } else {
    nlm.out <- runNLMChooseReg(X,y,grad=useGradient,useReg=useReg,lambda=lambda)
  }

  # set up the grid for visualizing the decision boundary
  thetaVec <- nlm.out$estimate
  
  col.names <-paste0("X",seq_along(names(X))-1)
  colnames(df) <- c(col.names,"y")
  grid.df <- data.frame(expand.grid(X1=seq(min(df$X1),max(df$X1),by=stepSize),X2=seq(min(df$X2),max(df$X2),by=stepSize)))
  if(!useEllipse && ellipseTransform){
   
    grid.df <- as.data.frame(cbind(rep(1,dim(grid.df)[1]), grid.mat))
    grid.df <- colsTransformEllipse(grid.df,k=3,mainEffects=FALSE)
    colnames(grid.df) <- col.names
  } else {
    grid.mat <- mapFeature(as.matrix(grid.df),degree=degree)
    grid.df <- as.data.frame(cbind(rep(1,dim(grid.df)[1]), grid.mat))
    colnames(grid.df) <- col.names
  }  
  
  
  # to get the probability of y=1 | X foreach example, we need to apply the learned theta
  # for each algorithm, this is done slightly differently...
  if(useEllipse){
    df$logLik <- apply(df[,1:3],1,function(row)t(row) %*% matrix(thetaVec,3) %*% row)
    stats.list <- testThetaStatsEllipse(X[,1:3], y, thetaVec)
    grid.df$logLik <- apply(grid.df[,1:3],1,function(row)t(row) %*% matrix(thetaVec,3) %*% row)
    grid.df$predict <- ifelse(grid.df$logLik < 0, 0, 1 )
  } else {
    df$logLik <- apply(X,1,function(row)t(row) %*% thetaVec)
    stats.list <- testThetaStats(df[col.names], y, thetaVec)
    grid.df$logLik <- apply(grid.df[col.names],1,function(row)t(row) %*% thetaVec)
    grid.df$predict <- ifelse(grid.df$logLik < 0, 0, 1 )
  }
  df$predict <- ifelse(df$logLik < 0, 0, 1 )
  # create grid for the prediction
  
  
  # figure out the title, with stats, dataset info, and lambda info(if regularized)
  stats.title <- with(stats.list,paste("\nAUC=",round(AUC,digits=2)," F1=",round(F1,digits=2),
                                       " sens=",round(sens,digits=2)," prec=",round(prec,digits=2),"\n",sep=""))
  title <- paste("mlclass ex 2 training dataset=",exprName," poly degree=",degree,stats.title,"\ncost fn = ",round(nlm.out$minimum,digits=2),sep="")
  title <- ifelse(useReg, paste(title, "\nreg'd w/ lambda =",lambda,sep=""),title)
  print(dim(grid.df)) 
  print(dim(df)) 
  ggplot(grid.df[c("X1","X2","predict")], aes(x= X1,y= X2,fill=factor(predict))) + 
    geom_raster(alpha=I(0.4)) +
    theme_bw() +
    layer( data = df,
           mapping = aes(x= X1,y= X2,color=factor( y),size=3),
           geom = "point" ) +
    ggtitle(title)
  
  ggsave(file=outputFile,height=8,width=8)
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
                               exprName = paste("ex2data",i,sep=""),useGradient=TRUE,useReg = FALSE)
    }
    
  }
  
  lambda.vals <- c(0,1,100)
  for(l in lambda.vals){
    dataset.file <- paste(getFullPath("data/logRegTest/mlclass-ex2/ex2data"), 2, ".txt",sep="")
    output.file <- paste(getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data2"),"_lambda=",l,"_degree=6",".pdf",sep="")
    plotPredictionML(inputFile=dataset.file,
                     outputFile=output.file,
                     degree=6,
                     exprName = paste("ex2data",2,sep=""),useGradient=TRUE,
                     useReg = TRUE, lambda = l, ySwap=FALSE, scaleData=FALSE,useEllipse=FALSE)
  }
  lambda.vals.ellipse <- c(0,1,10,50,100)
  # ml algo with poly degree 2
  for(l in lambda.vals.ellipse){
    dataset.file <- paste(getFullPath("data/logRegTest/mlclass-ex2/ex2data"), 2, ".txt",sep="")
    output.file <- paste(getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ex2data2"),"_lambda=",l,"_degree=2",".pdf",sep="")
    plotPredictionML(inputFile=dataset.file,
                     outputFile=output.file,
                     degree=2,
                     exprName = paste("ex2data",i,sep=""),useGradient=TRUE,
                     useReg = TRUE, lambda = l, ySwap=FALSE, scaleData=FALSE,useEllipse=FALSE,
                     msg = "ml algo with poly degree 2")
  }
  # elliptical algo 
  for(l in lambda.vals.ellipse){
    dataset.file <- paste(getFullPath("data/logRegTest/mlclass-ex2/ex2data"), 2, ".txt",sep="")
    output.file <- paste(getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/ellipse_ex2data2"),"_lambda=",l,"_degree=2",".pdf",sep="")
    plotPredictionML(inputFile=dataset.file,
                     outputFile=output.file,
                     degree=2,
                     exprName = paste("ex2data",i,sep=""),useGradient=TRUE,
                     useReg = TRUE, lambda = l, ySwap=FALSE, scaleData=FALSE,useEllipse=TRUE,
                     msg = "elliptical algo ")
  }
  # ml algo with elliptical function transform("colsTransformEllipse")
  for(l in lambda.vals.ellipse){
    dataset.file <- paste(getFullPath("data/logRegTest/mlclass-ex2/ex2data"), 2, ".txt",sep="")
    output.file <- paste(getFullPath("plots/fullAnalysisExperiment/test/logReg/mlclass/colsTransform_ex2data2"),"_lambda=",l,"_degree=2",".pdf",sep="")
    plotPredictionMLEllipseTransform(inputFile=dataset.file,
                     outputFile=output.file,
                     degree=2,
                     exprName = paste("ex2data",i,sep=""),useGradient=TRUE,
                     useReg = TRUE, lambda = l, ySwap=FALSE, scaleData=FALSE,useEllipse=FALSE,
                     #ellipseTransform=TRUE,
                     msg="ml-algo w/ elliptical trans")
  }
    
}

colsTransformEllipse <- function(X,k,mainEffects=FALSE){
  current <- k + 1
  for(i in 1:k){
    for (j in 1:k){
      X <- cbind(X,X[,i] * X[,j]) 
    }
  }
  if (mainEffects){
   X  
  }
  else{
    X[,(k+1):dim(X)[2]]
  }
  
}




