rm(list=unlist(ls()))

library(Rcpp)
library(inline)
library(RcppArmadillo)
# p 26 

src1 <- '
  Rcpp::NumericVector xa(a);
  Rcpp::NumericVector xb(b);
  int n_xa = xa.size(), n_xb = xb.size();

Rcpp::NumericVector xab(n_xa + n_xb + 1);
for (int i = 0; i < n_xa; i++)
  for (int j = 0; j < n_xb; j++)
    xab[i + j] += xa[i] * xb[j];
  return xab;
'






fun1 = cxxfunction(signature(a="numeric", b="numeric"), src1, plugin="Rcpp")
fun1(1:4,2:5)



# armadillo example:
# http://arma.sourceforge.net/docs.html

src2 <- '
  Rcpp::NumericMatrix Xr(Xs);
  Rcpp::NumericVector yr(ys);
  int n = Xr.nrow(),  k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false);
  arma::colvec y(yr.begin(), yr.size(), false);
  arma::mat sim  = 1 / (trans(X)*X);
  arma::mat eigenVectors;
  arma::colvec eigenValues;
  eig_sym(eigenValues, eigenVectors, sim);
  return Rcpp::List::create(Rcpp::Named("vector") = eigenValues );
'

f2 = cxxfunction(signature(Xs="numeric", ys="numeric"), src2, plugin="RcppArmadillo")

#k = 10000;f2(Xs=matrix(runif(k*k),k,k),runif(k))


src4 <- '
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
f4 = cxxfunction(signature(Xs="numeric", ys="numeric"), src4, plugin="RcppArmadillo")

a= f4(Xs = matrix(runif(500*40),500), ys = runif(500))



X = as.matrix(iris[,1:4])
y = runif(dim(X)[1])
ic = f4(Xs = X, ys = y)
yp = as.matrix(ic$y,1)


X = as.matrix(iris[,1:4])
S = 1/ (X %*% t(X))
e = eigen(S)
ye = e$vector




src3 = '
  int n = 100;
  int k = 30;
  arma::mat X = arma::randu<arma::mat>(100,30);
//  X.randu(n,k);
  arma::mat sim = (X*trans(X));
  sim = arma::pow(sim,-1);
  arma::mat eigenVectors;
  arma::colvec eigenValues;
  eig_sym(eigenValues, eigenVectors, sim, "dc");
  //eigenValues.max(index);
  arma::colvec output = eigenVectors.row(1);
  return Rcpp::List::create(Rcpp::Named("eigenVector") = output );
';f3 = cxxfunction(signature(), src3, plugin="RcppArmadillo");f3()








