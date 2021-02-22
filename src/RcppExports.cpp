// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_ind
arma::vec get_ind(arma::mat b_mat, int n, int H);
RcppExport SEXP _pjnormarcpp_get_ind(SEXP b_matSEXP, SEXP nSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type b_mat(b_matSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ind(b_mat, n, H));
    return rcpp_result_gen;
END_RCPP
}
// generate_b_cube
arma::cube generate_b_cube(arma::vec ind, int k, int H, int n);
RcppExport SEXP _pjnormarcpp_generate_b_cube(SEXP indSEXP, SEXP kSEXP, SEXP HSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_b_cube(ind, k, H, n));
    return rcpp_result_gen;
END_RCPP
}
// get_sigmamat_full
arma::mat get_sigmamat_full(arma::mat sigmamat_conditional, arma::vec gamma, int k);
RcppExport SEXP _pjnormarcpp_get_sigmamat_full(SEXP sigmamat_conditionalSEXP, SEXP gammaSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigmamat_conditional(sigmamat_conditionalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sigmamat_full(sigmamat_conditional, gamma, k));
    return rcpp_result_gen;
END_RCPP
}
// subtract_y
arma::mat subtract_y(arma::cube phi_arr, arma::vec ind, arma::mat y, int n);
RcppExport SEXP _pjnormarcpp_subtract_y(SEXP phi_arrSEXP, SEXP indSEXP, SEXP ySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type phi_arr(phi_arrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(subtract_y(phi_arr, ind, y, n));
    return rcpp_result_gen;
END_RCPP
}
// subtract_b_cube
arma::cube subtract_b_cube(arma::cube b_cube, arma::cube phi_arr, arma::vec ind, int n);
RcppExport SEXP _pjnormarcpp_subtract_b_cube(SEXP b_cubeSEXP, SEXP phi_arrSEXP, SEXP indSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type b_cube(b_cubeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type phi_arr(phi_arrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(subtract_b_cube(b_cube, phi_arr, ind, n));
    return rcpp_result_gen;
END_RCPP
}
// generatemu_posterior
List generatemu_posterior(arma::mat mu_prior_precision, arma::cube precision_full_arr, arma::vec ind, arma::vec mu_prior_mu, arma::cube subtracted_b_cube, arma::mat subtracted_y, int n);
RcppExport SEXP _pjnormarcpp_generatemu_posterior(SEXP mu_prior_precisionSEXP, SEXP precision_full_arrSEXP, SEXP indSEXP, SEXP mu_prior_muSEXP, SEXP subtracted_b_cubeSEXP, SEXP subtracted_ySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mu_prior_precision(mu_prior_precisionSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type precision_full_arr(precision_full_arrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_prior_mu(mu_prior_muSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type subtracted_b_cube(subtracted_b_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type subtracted_y(subtracted_ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(generatemu_posterior(mu_prior_precision, precision_full_arr, ind, mu_prior_mu, subtracted_b_cube, subtracted_y, n));
    return rcpp_result_gen;
END_RCPP
}
// center_y
arma::mat center_y(arma::mat subtracted_y, arma::cube subtracted_b_cube, arma::vec mu, int n);
RcppExport SEXP _pjnormarcpp_center_y(SEXP subtracted_ySEXP, SEXP subtracted_b_cubeSEXP, SEXP muSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type subtracted_y(subtracted_ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type subtracted_b_cube(subtracted_b_cubeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(center_y(subtracted_y, subtracted_b_cube, mu, n));
    return rcpp_result_gen;
END_RCPP
}
// generategamma_posterior
List generategamma_posterior(arma::vec ind, arma::cube sigmamat_conditional_arr, arma::mat gamma_prior_precision, arma::vec gamma_prior_mu, arma::mat centered_y, int h, int k);
RcppExport SEXP _pjnormarcpp_generategamma_posterior(SEXP indSEXP, SEXP sigmamat_conditional_arrSEXP, SEXP gamma_prior_precisionSEXP, SEXP gamma_prior_muSEXP, SEXP centered_ySEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigmamat_conditional_arr(sigmamat_conditional_arrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma_prior_precision(gamma_prior_precisionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma_prior_mu(gamma_prior_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type centered_y(centered_ySEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(generategamma_posterior(ind, sigmamat_conditional_arr, gamma_prior_precision, gamma_prior_mu, centered_y, h, k));
    return rcpp_result_gen;
END_RCPP
}
// diff_centered_y
arma::mat diff_centered_y(arma::mat centered_y, arma::mat gamma_mat, arma::vec ind, int k, int n);
RcppExport SEXP _pjnormarcpp_diff_centered_y(SEXP centered_ySEXP, SEXP gamma_matSEXP, SEXP indSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type centered_y(centered_ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gamma_mat(gamma_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(diff_centered_y(centered_y, gamma_mat, ind, k, n));
    return rcpp_result_gen;
END_RCPP
}
// get_w
arma::mat get_w(arma::mat y, arma::vec mu, arma::vec ind, int H, int k, int n);
RcppExport SEXP _pjnormarcpp_get_w(SEXP ySEXP, SEXP muSEXP, SEXP indSEXP, SEXP HSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type H(HSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(get_w(y, mu, ind, H, k, n));
    return rcpp_result_gen;
END_RCPP
}
// generatephi_posterior
List generatephi_posterior(arma::vec ind1, arma::mat w, arma::cube precision_full_arr, arma::cube phi_arr, arma::mat phi_prior_precision, arma::vec phi_prior_mu, int k, int h, int kk);
RcppExport SEXP _pjnormarcpp_generatephi_posterior(SEXP ind1SEXP, SEXP wSEXP, SEXP precision_full_arrSEXP, SEXP phi_arrSEXP, SEXP phi_prior_precisionSEXP, SEXP phi_prior_muSEXP, SEXP kSEXP, SEXP hSEXP, SEXP kkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type precision_full_arr(precision_full_arrSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type phi_arr(phi_arrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_prior_precision(phi_prior_precisionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi_prior_mu(phi_prior_muSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type kk(kkSEXP);
    rcpp_result_gen = Rcpp::wrap(generatephi_posterior(ind1, w, precision_full_arr, phi_arr, phi_prior_precision, phi_prior_mu, k, h, kk));
    return rcpp_result_gen;
END_RCPP
}
// gen_r
double gen_r(arma::vec ind, arma::mat mu_present_mat, arma::cube precision_full_arr, arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn, int n, double r, int k);
RcppExport SEXP _pjnormarcpp_gen_r(SEXP indSEXP, SEXP mu_present_matSEXP, SEXP precision_full_arrSEXP, SEXP phi_arrSEXP, SEXP sphericaldataSEXP, SEXP augmenteddata_rawSEXP, SEXP nnSEXP, SEXP nSEXP, SEXP rSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_present_mat(mu_present_matSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type precision_full_arr(precision_full_arrSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type phi_arr(phi_arrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sphericaldata(sphericaldataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type augmenteddata_raw(augmenteddata_rawSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_r(ind, mu_present_mat, precision_full_arr, phi_arr, sphericaldata, augmenteddata_raw, nn, n, r, k));
    return rcpp_result_gen;
END_RCPP
}
// gen_r_group
double gen_r_group(arma::vec ind, arma::vec ind1, arma::mat mu_present_mat, arma::cube precision_full_arr, arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn, int n, double r, int k);
RcppExport SEXP _pjnormarcpp_gen_r_group(SEXP indSEXP, SEXP ind1SEXP, SEXP mu_present_matSEXP, SEXP precision_full_arrSEXP, SEXP phi_arrSEXP, SEXP sphericaldataSEXP, SEXP augmenteddata_rawSEXP, SEXP nnSEXP, SEXP nSEXP, SEXP rSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind1(ind1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_present_mat(mu_present_matSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type precision_full_arr(precision_full_arrSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type phi_arr(phi_arrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sphericaldata(sphericaldataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type augmenteddata_raw(augmenteddata_rawSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_r_group(ind, ind1, mu_present_mat, precision_full_arr, phi_arr, sphericaldata, augmenteddata_raw, nn, n, r, k));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_chain_pjnormar1
List mcmc_chain_pjnormar1(arma::mat sphericaldata, arma::mat b_mat, int iterations);
RcppExport SEXP _pjnormarcpp_mcmc_chain_pjnormar1(SEXP sphericaldataSEXP, SEXP b_matSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sphericaldata(sphericaldataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b_mat(b_matSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_chain_pjnormar1(sphericaldata, b_mat, iterations));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_chain_pjnormar1_group
List mcmc_chain_pjnormar1_group(arma::mat sphericaldata, arma::uvec protein_first, arma::mat b_mat, int iterations);
RcppExport SEXP _pjnormarcpp_mcmc_chain_pjnormar1_group(SEXP sphericaldataSEXP, SEXP protein_firstSEXP, SEXP b_matSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sphericaldata(sphericaldataSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type protein_first(protein_firstSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b_mat(b_matSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_chain_pjnormar1_group(sphericaldata, protein_first, b_mat, iterations));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::vec mvrnormArma(int k, arma::vec mu, arma::mat sigma);
RcppExport SEXP _pjnormarcpp_mvrnormArma(SEXP kSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(k, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// wishartArma
arma::mat wishartArma(int k, int n, arma::mat sigma);
RcppExport SEXP _pjnormarcpp_wishartArma(SEXP kSEXP, SEXP nSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(wishartArma(k, n, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pjnormarcpp_get_ind", (DL_FUNC) &_pjnormarcpp_get_ind, 3},
    {"_pjnormarcpp_generate_b_cube", (DL_FUNC) &_pjnormarcpp_generate_b_cube, 4},
    {"_pjnormarcpp_get_sigmamat_full", (DL_FUNC) &_pjnormarcpp_get_sigmamat_full, 3},
    {"_pjnormarcpp_subtract_y", (DL_FUNC) &_pjnormarcpp_subtract_y, 4},
    {"_pjnormarcpp_subtract_b_cube", (DL_FUNC) &_pjnormarcpp_subtract_b_cube, 4},
    {"_pjnormarcpp_generatemu_posterior", (DL_FUNC) &_pjnormarcpp_generatemu_posterior, 7},
    {"_pjnormarcpp_center_y", (DL_FUNC) &_pjnormarcpp_center_y, 4},
    {"_pjnormarcpp_generategamma_posterior", (DL_FUNC) &_pjnormarcpp_generategamma_posterior, 7},
    {"_pjnormarcpp_diff_centered_y", (DL_FUNC) &_pjnormarcpp_diff_centered_y, 5},
    {"_pjnormarcpp_get_w", (DL_FUNC) &_pjnormarcpp_get_w, 6},
    {"_pjnormarcpp_generatephi_posterior", (DL_FUNC) &_pjnormarcpp_generatephi_posterior, 9},
    {"_pjnormarcpp_gen_r", (DL_FUNC) &_pjnormarcpp_gen_r, 10},
    {"_pjnormarcpp_gen_r_group", (DL_FUNC) &_pjnormarcpp_gen_r_group, 11},
    {"_pjnormarcpp_mcmc_chain_pjnormar1", (DL_FUNC) &_pjnormarcpp_mcmc_chain_pjnormar1, 3},
    {"_pjnormarcpp_mcmc_chain_pjnormar1_group", (DL_FUNC) &_pjnormarcpp_mcmc_chain_pjnormar1_group, 4},
    {"_pjnormarcpp_mvrnormArma", (DL_FUNC) &_pjnormarcpp_mvrnormArma, 3},
    {"_pjnormarcpp_wishartArma", (DL_FUNC) &_pjnormarcpp_wishartArma, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_pjnormarcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}