#include <armadillo>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;


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
/*





f = cxxfunction(signature(xs="numeric",ys="numeric"),
      file = paste( readLines(getFullPath("/analysis/rcpp_src/arm_test.cpp" ), collapse = "\n" ),
      plugin="Rcpp"))











