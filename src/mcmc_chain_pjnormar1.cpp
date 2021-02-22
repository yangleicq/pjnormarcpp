#include <RcppArmadillo.h>
using namespace Rcpp;
#include "innerFunctions.h"
#include "random.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List mcmc_chain_pjnormar1(arma::mat sphericaldata, arma::mat b_mat, int iterations=5000){
  //dimension
  int k = sphericaldata.n_cols;
  int kk;
  //regime number
  int H = b_mat.n_cols;
  int h;
  //sample size
  int n = sphericaldata.n_rows;
  int nn;
  //right now only consider AR(1)
  int p = 1;


  //values used for simulation

  arma::vec n_vec = sum(b_mat, 0).t();

  arma::vec ind = get_ind(b_mat, n, H);
  arma::vec ind1 = ind;
  ind1(0) = 999;

  arma::cube b_cube = generate_b_cube(ind, k, H, n);

  //Initial values
  // mu is mean vector of length (k*H)
  arma::vec mu_present(k*H); mu_present.zeros();
  arma::mat mu_present_mat(k,H); mu_present_mat.zeros();

  // gamma_mat is a (k-1) by H matrix with [,h] as gamma_present for regime h
  arma::mat gamma_mat((k-1),H); gamma_mat.zeros();
  // sigma_conditional_arr is (k-1) by (k-1) by H array with [,,h]
  // as sigma_conditional_present for regime h
  arma::cube sigmamat_conditional_arr((k-1),(k-1),H);
  arma::mat diagk_1((k-1),(k-1)); diagk_1.eye();
    for(int h=0; h<H; h++)
        sigmamat_conditional_arr.slice(h) = diagk_1;
  //compute the precision_full_arr from sigmamat_conditional_arr and gamma_mat
  arma::cube precision_full_arr(k,k,H);
  arma::mat sigmamat_full;
  for(h=0; h<H; h++){
    sigmamat_full = get_sigmamat_full(sigmamat_conditional_arr.slice(h), gamma_mat.col(h),k);
    precision_full_arr.slice(h) = inv_sympd(sigmamat_full);
    }

  // phi_arr is k by k by H array with [,,h]
  // as phi_present for regime h
  arma::cube phi_arr(k,k,H); phi_arr.zeros();

  // vector lengths
  arma::vec r_present(n);  r_present.ones();

  arma::mat augmenteddata_raw(n,k);
  for(kk=0; kk<k; kk++)
    augmenteddata_raw.col(kk)=sphericaldata.col(kk)%r_present;

  //declaration for variables used in simulation
  arma::cube subtracted_b_cube;
  arma::mat subtracted_y;
  List mu_posterior;
  arma::mat centered_y;
  List gamma_posterior;
  arma::mat difference;
  arma::uvec seq_vech;
  arma::mat sigmamat_contitional_posterior_mat;
  arma::mat sigmamat_conditional_inverse;
  arma::mat w;
  List phi_posterior;


  //containers for sampling variables
  arma::mat mu_all(iterations, k*H);
  arma::cube gamma_mat_all((k-1), H, iterations);
  List sigmamat_conditional_arr_all(iterations);
  List phi_arr_all(iterations);
  arma::cube augmenteddata_raw_all(n,k,iterations);
  arma::mat r_all(iterations, n);
  //priors
  //we use the same prior for all H regimes for simplicity
  //mu~Normal
  arma::vec mu_prior_mu(k*H); mu_prior_mu.zeros();
  arma::mat mu_prior_sigmamat(k*H,k*H);mu_prior_sigmamat.eye();
  mu_prior_sigmamat = mu_prior_sigmamat*1e5;
  arma::mat mu_prior_precision = inv_sympd(mu_prior_sigmamat);
  //gamma~Normal
  arma::vec gamma_prior_mu(k-1); gamma_prior_mu.zeros();
  arma::mat gamma_prior_sigmamat(k-1,k-1); gamma_prior_sigmamat.eye();
  gamma_prior_sigmamat = gamma_prior_sigmamat*1e5;
  arma::mat gamma_prior_precision = inv_sympd(gamma_prior_sigmamat);
  //sigma~Wishart
  arma::vec sigmamat_conditional_prior_df(H); sigmamat_conditional_prior_df.ones();
  sigmamat_conditional_prior_df = sigmamat_conditional_prior_df*4;
  arma::mat sigmamat_conditional_prior_mat(k-1,k-1); sigmamat_conditional_prior_mat.eye();
  //phi~Normal
  arma::vec phi_prior_mu(k); phi_prior_mu.zeros();
  arma::mat phi_prior_sigmamat(k,k); phi_prior_sigmamat.eye();
  phi_prior_sigmamat = phi_prior_sigmamat*1e5;
  arma::mat phi_prior_precision = inv_sympd(phi_prior_sigmamat);

  //set seed
//  arma::arma_rng::set_seed(seed);


  for(int iter=0; iter<iterations; iter++)
  {
    //sampling mu
    subtracted_b_cube = subtract_b_cube(b_cube, phi_arr, ind, n);
    subtracted_y = subtract_y(phi_arr, ind, augmenteddata_raw, n);
    mu_posterior = generatemu_posterior(mu_prior_precision, precision_full_arr,
                                        ind, mu_prior_mu, subtracted_b_cube, subtracted_y, n);
    mu_present = mvrnormArma(k*H,mu_posterior["mu_posterior_mu"], mu_posterior["mu_posterior_sigmamat"]);
    for(h=0; h<H; h++){
      mu_present_mat.col(h) = mu_present(arma::span(h*k,h*k+k-1));
    }
    mu_all.row(iter) = mu_present.t();

    //sampling gamma_h

    centered_y = center_y(subtracted_y,subtracted_b_cube, mu_present, n);
    for(h=0; h<H; h++){
      gamma_posterior = generategamma_posterior(ind, sigmamat_conditional_arr,gamma_prior_precision, gamma_prior_mu, centered_y, h, k);
      gamma_mat.col(h) = mvrnormArma((k-1), gamma_posterior["gamma_posterior_mu"], gamma_posterior["gamma_posterior_sigmamat"]);
      }

    gamma_mat_all.slice(iter) = gamma_mat;

    //sampling sigmamat_conditional

    difference = diff_centered_y(centered_y, gamma_mat, ind, k, n);
    for(h=0; h<H; h++){
      seq_vech = find(ind==h);
      sigmamat_contitional_posterior_mat = sigmamat_conditional_prior_mat+difference.rows(seq_vech).t()*difference.rows(seq_vech);
      sigmamat_conditional_inverse = wishartArma((k-1), n_vec(h)+sigmamat_conditional_prior_df(h), inv_sympd(sigmamat_contitional_posterior_mat));
      sigmamat_conditional_arr.slice(h) = inv_sympd(sigmamat_conditional_inverse);
    }
      sigmamat_conditional_arr_all[iter] = sigmamat_conditional_arr;

    //compute sigmamat_full
    for(h=0; h<H; h++){
      sigmamat_full = get_sigmamat_full(sigmamat_conditional_arr.slice(h), gamma_mat.col(h), k);
      precision_full_arr.slice(h) = inv_sympd(sigmamat_full);
    }

    //sampling for phi
    w = get_w(augmenteddata_raw, mu_present, ind, H, k, n);

    for(h=0; h<H; h++)
      for(kk=0; kk<k; kk++){
        phi_posterior = generatephi_posterior(ind1, w, precision_full_arr, phi_arr, phi_prior_precision, phi_prior_mu,k,h,kk);
        phi_arr.slice(h).col(kk) = mvrnormArma(k, phi_posterior["phi_posterior_mu"], phi_posterior["phi_posterior_sigmamat"]);
      }
      phi_arr_all[iter] = phi_arr;

    //sampling for r
    for(nn=0; nn<n; nn++){
     r_present(nn) = gen_r(ind, mu_present_mat, precision_full_arr,
            phi_arr, sphericaldata, augmenteddata_raw, nn, n, r_present(nn), k);
    augmenteddata_raw.row(nn) = r_present(nn)*sphericaldata.row(nn);
    }
    r_all.row(iter) = r_present.t();
//    std::cout << iter << std::endl;
  }
  return List::create(Rcpp::Named("mu_all") = mu_all,
                      Rcpp::Named("gamma_mat_all") = gamma_mat_all,
                      Rcpp::Named("sigmamat_conditional_arr_all") = sigmamat_conditional_arr_all,
                      Rcpp::Named("phi_arr_all") = phi_arr_all,
                      Rcpp::Named("r_all") = r_all);

}

