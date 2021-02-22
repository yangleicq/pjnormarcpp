#ifndef INNERFUNCTIONS_H
#define INNERFUNCTIONS_H

#include <RcppArmadillo.h>
using namespace Rcpp;
arma::vec get_ind(arma::mat b_mat, int n, int H);

arma::cube generate_b_cube(arma::vec ind, int k, int H, int n);

arma::mat get_sigmamat_full(arma::mat sigmamat_conditional, arma::vec gamma, int k);

arma::mat subtract_y(arma::cube phi_arr,arma::vec ind, arma::mat y, int n);

arma::cube subtract_b_cube(arma::cube b_cube, arma::cube phi_arr, arma::vec ind, int n);

List generatemu_posterior(arma::mat mu_prior_precision, arma::cube precision_full_arr,
                          arma::vec ind, arma::vec mu_prior_mu, arma::cube subtracted_b_cube, arma::mat subtracted_y, int n);

arma::mat center_y(arma::mat subtracted_y, arma::cube subtracted_b_cube, arma::vec mu, int n);

List generategamma_posterior(arma::vec ind, arma::cube sigmamat_conditional_arr,
                               arma::mat gamma_prior_precision, arma::vec gamma_prior_mu, arma::mat centered_y, int h, int k);

arma::mat diff_centered_y(arma::mat centered_y, arma::mat gamma_mat, arma::vec ind, int k, int n);

arma::mat get_w(arma::mat y, arma::vec mu, arma::vec ind, int H, int k, int n);

List generatephi_posterior(arma::vec ind1, arma::mat w, arma::cube precision_full_arr, arma::cube phi_arr,
                             arma::mat phi_prior_precision, arma::vec phi_prior_mu, int k, int h, int kk);

double gen_r(arma::vec ind, arma::mat mu_present_mat, arma::cube precision_full_arr,
             arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn,
             int n, double r, int k);

double gen_r_group(arma::vec ind, arma::vec ind1, arma::mat mu_present_mat, arma::cube precision_full_arr,
                   arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn,
                   int n, double r, int k);

#endif
