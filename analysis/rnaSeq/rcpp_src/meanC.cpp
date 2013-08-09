#include <Rcpp.h>
#include <armadillo>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i] / n;
  }
  return total;
}

// [[Rcpp::export]]
SEXP testM(SEXP ys, SEXP xs){
  // testM :: fn rVector -> rMatrix -> rVector
  //returns a matrix of double
try{
    // get Rcpp objects from the two S expression
    Rcpp::NumericVector yr(ys);
    Rcpp::NumericMatrix Xr(xs);
    
    int n = Xr.nrow(); int k = Xr.ncol();
    
    //pass the Rcpp datastructs into armadillo
    arma::mat X(Xr.begin(), n, k, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec out = trans(X) * y;
    
    return wrap(out);
  } catch( std::exception &ex){
    forward_exception_to_r( ex );
  } catch(...){
    ::Rf_error( "C++ exception (unknown reason) ~ arm_test.cpp\n");
  }
return R_NilValue;   
}


/**************************************************
 * 
 * Demonstrantion of using another function defined in the namespace...
 * 
 */



int adder_internal(int x, int y){
  int out = (x + y)/ (x * y);
  return out;
}

// [[Rcpp::export]]
SEXP applyM(SEXP vec, int x, int y){
  // testM :: fn rVector -> rMatrix -> rVector
  //returns a vector with two ints added to it (map #(+ % x y) vec)
  try{
    // get Rcpp objects from the two S expression
    Rcpp::NumericVector yr(vec);
    int n = yr.size();
    
    NumericVector out(n);

    for (int i = 0; i < n; ++i) {
      out[i] = out[i] + x + y;
    }

    return out;

    int j = adder_internal( x , y);
    return wrap(yr);
  } catch( std::exception &ex){
    forward_exception_to_r( ex );
  } catch(...){
    ::Rf_error( "C++ exception (unknown reason) ~ arm_test.cpp\n");
  }
return R_NilValue;   
}



/*
// [[Rcpp::export]]
extern "C" SEXP testM(SEXP ys, SEXP xs){
  try{
    // get Rcpp objects from the two S expression
    Rcpp::NumericVector yr(ys);
    Rcpp::NumericMatrix Xr(xs);
    
    // get the dimensions of Xr
    int n = Xr.nrow(); int k = Xr.ncol();
    
    //pass the Rcpp datastructs into armadillo
    arma::mat X(Xr.begin(), n, k, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    arma::colvec out = trans(X) * y;
    return Rcpp::List::create(Rcpp::Named("vector") = out );
    
    //catch some expections here to pass back to R...
  } catch( std::exception &ex){
    forward_exception_to_r( ex );
  } catch(...){
    ::Rf_error( "C++ exception (unknown reason) ~ arm_test.cpp\n");
  }
  return R_NilValue;  
  
}
*/
