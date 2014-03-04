
# http://leizhengdeng.blogspot.com/2012/07/order-statistics-for-combining-various.html

Q <- function(r){
  r <- sort(r)
  N <- length(r)
  if (N==1){
    return (r)
  }
  s <- 0
  for (i in 1:N){
    del.element.idx <- N-i+1
    if (N==i){
      s = s+(r[N-i+1])*Q(r[-del.element.idx])
    }else{
      s = s+(r[N-i+1]-r[N-i])*Q(r[-del.element.idx])
    }
  }
  s <- factorial(N)*s
  return (s)
}

#r <- c(0.2, 0.05, 0.2)
#Q(r)
appendT <- function(tabs){cat(rep('  ',tabs))} 

corByRow <- function(df1){
A <- as.matrix(df1)
B <- as.matrix(df1)
sapply(seq.int(dim(A)[1]), function(i) cor(A[i,], B[i,]))
}


Q.test <- function(r,t){
  appendT(t);cat("Calling Q (top), r=",r,"\n")
  r <- sort(r)
  N <- length(r)
  if (N==1){
    appendT(t);cat("returning N==1 ->",r,"\n")
    return (r)
  }
  s <- 0
  for (i in 1:N){
    del.element.idx <- N-i+1
    appendT(t);cat("del index -> ",del.element.idx,"\n")
    t.next <- t + 1
    if (N==i){
      s = s+(r[N-i+1])*Q.test(r[-del.element.idx],t.next)
    }else{
      s = s+(r[N-i+1]-r[N-i])*Q.test(r[-del.element.idx],t.next)
    }
  }
  #s <- factorial(N)*s
  return (s)
}


d<-matrix(c(sample(1:10),sample(1:10),sample(1:10),sample(1:10)), nrow=10, byrow=F)
rownames(d) <- c(paste("gene", seq(1:10), sep=""))
colnames(d) <- c("gene.exp", "methyl", "CNV", "Mutation")
d
r = d/10
ranking <- apply(r, 1, Q)
sort(ranking)

#http://bioops.info/2013/08/p-value-order-statistics/
nsquare = 200
df<- data.frame(trial=seq_len(nsquare*nsquare))
df$cond.2 <- c(cor(matrix(runif(2*nsquare),nrow=2)))
df$cond.5 <- c(cor(matrix(runif(5*nsquare),nrow=5)))
df$cond.15 <- c(cor(matrix(runif(15*nsquare),nrow=15)))
df$cond.50 <- c(cor(matrix(runif(50*nsquare),nrow=50)))
df.melt <- melt(df,id.var="trial")
colnames(df.melt) <- c("trial", "conditions", "R")
ggplot(df.melt, aes(x=R,fill=conditions))+geom_density(alpha=I(0.4)) + theme_bw()


