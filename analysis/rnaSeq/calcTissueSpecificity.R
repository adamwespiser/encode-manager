
# Created by Adam Wespiser


# entropy = - sum ( x * log(x))
# if x[n] == 0, then log(x[n]) will be NA, we will set the contribution of this point to the total entropy at 0
entropyVec <- function(x){
  x <- x / sum(x)
  y <- log2(x) * x
  y[is.na(y)] <- 0
  sum(y) * -1
}


JS = function(x,y){
   x <- x / sum(x)
   y <- y / sum(y)
  
  a <- (x + y)/2
  z <- entropyVec(a) - ((entropyVec(x)-entropyVec(y))/2)
  z
}

JSsp <- function(e1,e2){
  1 - sqrt(JS(e1,e2))
}


# Use this function to calculate cell-type specificity of an expression vector
calcTissSpec <- function(x,cols=seq_along(x)){
  x <- sapply(x,as.numeric)
  profile <- diag(length(cols))
  specEn <- sapply(1:length(cols),function(y)JSsp(x,profile[y,]))
  
  max(specEn)
}