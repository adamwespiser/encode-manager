library(Rcpp)
library(inline)
library(RcppArmadillo)




srcCalcEigenCpp <- '
Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericVector yr(ys);
int m = as<int>(mat_flag);
int n = Xr.nrow(),  k = Xr.ncol();
arma::mat X(Xr.begin(), n, k, false);
arma::colvec y(yr.begin(), yr.size(), false);
arma::colvec y_prev = y;
arma::colvec y_temp = y;

arma::mat sim =  X*trans(X);
double found_zero = 0;

if ( m == 0){
sim = arma::pow(sim,-1);
m = -2;
}
else {

m = -1;
int i;

for(i = 0; i < n; i++){
sim(i,i) = 1/sim(i,i);
int j;
for (j = i + 1; j < n; j++){

m = j *10000 + i;
double val = sim(i,j);
if (abs(val) < 0.001){
//sim(i,j) = sim(j,i) = 0;	
} 
else {
sim(i,j) = sim(j,i) = pow(val,-1);
}
}
}
}

int i = 0;
double diff = 10000000.0;
double y_norm = 0;
while(diff > 1e-15 && i < 2000){
i++; 
y_prev = y;

y_temp = sim * y;
y_norm = norm(y_temp,2);
y = y_temp / y_norm ;
//diff = sum(abs(y_prev - y)) ;
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
Rcpp::Named("mat")       = m,
Rcpp::Named("foundZero")       = found_zero,
Rcpp::Named("diff")       = diff);

' 
# end of src4....
# ... 
CalcEigenCpp1 = cxxfunction(signature(Xs="numeric", ys="numeric",mat_flag="int"), srcCalcEigenCpp, plugin="RcppArmadillo")







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
     distance(i,j) = distance(j,i) = distValue;

  }
}

i = 0;
double diff = 10000000.0;
double y_norm = 0;
while(diff > 1e-15 && i < 2000){
i++; 
y_prev = y;

y_temp = sim * y;
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
cp = CalcEigenCppD(Xs=m,y=rep(1,dim(m)[1]))
rp = eigen(dist(m,upper=TRUE,diag=TRUE))
dm = dist(m,upper=TRUE,diag=TRUE)
as.vector(cp$y)
rp$vectors[,1]

if(abs(cov(cp$y,rp$vectors[,1])) > 0.95) {
  print("Vectors look to be the same")
  
}


