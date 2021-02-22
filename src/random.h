#ifndef RANDOM_H
#define RANDOM_H

#include <RcppArmadillo.h>
using namespace Rcpp;
#include <math.h>

arma::vec mvrnormArma(int k, arma::vec mu, arma::mat sigma);

arma::mat wishartArma(int k, int n, arma::mat sigma);



#endif
