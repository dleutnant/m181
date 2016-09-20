// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export(.run_constancy)]]
LogicalVector run_constancy(arma::vec x, int n, double delta_max, int method){
  
  int sz = x.size();
  
  LogicalVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
    if (method == 1) {
      // method 1 uses the original vector
      res[i+n-1] = std::abs(x.subvec(i, i + n - 1).max() - 
                            x.subvec(i, i + n - 1).min()) < delta_max;
      
    }

    if (method == 2) {
      // method 2 uses the diff vector
      res[i+n-1] = arma::all(x.subvec(i, i + n - 1) < delta_max);;
      
    }
    
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.run_mean)]]
NumericVector run_mean(arma::vec x, int n){
  
  int sz = x.size();
  
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
      res[i+n-1] = arma::mean(x.subvec(i, i + n - 1));
      
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.run_median)]]
NumericVector run_median(arma::vec x, int n){
  
  int sz = x.size();
  
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
    res[i+n-1] = arma::median(x.subvec(i, i + n - 1));
    
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.run_sd)]]
NumericVector run_sd(arma::vec x, int n){
  
  int sz = x.size();
  
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
    res[i+n-1] = arma::stddev(x.subvec(i, i + n - 1));
    
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.run_max)]]
NumericVector run_max(arma::vec x, int n){
  
  int sz = x.size();
  
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
    res[i+n-1] = x.subvec(i, i + n - 1).max();
    
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.run_min)]]
NumericVector run_min(arma::vec x, int n){
  
  int sz = x.size();
  
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
    res[i+n-1] = x.subvec(i, i + n - 1).min();
    
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.run_mad)]]
NumericVector run_mad(arma::vec x, int n, double scale_factor = 1.4826) {
  
  int sz = x.size();
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n+1); i++){
    
    res[i+n-1] = arma::median( abs(x.subvec(i, i + n - 1) - 
                               arma::median(x.subvec(i, i + n - 1))) ) * scale_factor;
    
  }
  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  return res;
}

// [[Rcpp::export(.count_sign_change)]]
NumericVector count_sign_change(arma::vec x, int n) {
  
  int sz = x.size();
  // diff to get gradients
  // sign to get sign of gradients
  // diff to get change of signs, i.e. everything but 0 is a sign change
  arma::vec x_new = diff(sign(diff(x)));
  
  NumericVector res(sz);
  
  for(int i = 0; i < (sz-n - 1); i++){
    // subset
    arma::vec x_new_sub = x_new.subvec(i, i + n - 1);
    // get all zeros
    arma::vec x_new_sub2 = x_new_sub.elem(arma::find(x_new_sub != 0 ));
    res[i+n-1] = x_new_sub2.n_elem;
  }

  // pad the first n-1 elements with NA
  std::fill(res.begin(), res.end()-sz+n-1, NA_REAL);
  
  // do not forget to add two NA's at the front and to remove the last two ele!
  
  return res;
}