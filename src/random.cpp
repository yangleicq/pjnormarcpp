#include <RcppArmadillo.h>
using namespace Rcpp;
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec mvrnormArma(int k, arma::vec mu, arma::mat sigma){
  arma::vec Y = arma::randn(k);
  return mu + arma::chol(sigma).t()*Y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat wishartArma(int k, int n, arma::mat sigma){
  arma::mat Y = arma::randn(k, n);
  return arma::chol(sigma).t()*Y*Y.t()*arma::chol(sigma);
}

