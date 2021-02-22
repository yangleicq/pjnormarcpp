#include <RcppArmadillo.h>
using namespace Rcpp;
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec get_ind(arma::mat b_mat, int n, int H) {
  arma::vec ind(n); ind.zeros();
  arma::rowvec seq_H = arma::linspace<arma::rowvec>(0,(H-1),H);
  for(int nn=0; nn<n; nn++){
    ind(nn)=sum(b_mat.row(nn)%seq_H);
  }
  return ind;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube generate_b_cube(arma::vec ind, int k, int H, int n){
  arma::cube b_cube(k,k*H,n); b_cube.zeros();
  arma::mat diagk(k,k); diagk.eye();
  for (int nn=0; nn<n; nn++){
    b_cube.slice(nn).submat(0,ind(nn)*k,k-1, (ind(nn)*k+k-1)) = diagk;
  }
  return b_cube;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_sigmamat_full(arma::mat sigmamat_conditional, arma::vec gamma, int k){
  arma::mat sigmamat_full(k,k); sigmamat_full.zeros();
  sigmamat_full.submat(0,0,(k-2),(k-2)) = sigmamat_conditional + gamma*gamma.t();
  sigmamat_full.submat(0,(k-1),(k-2),(k-1)) = gamma;
  sigmamat_full.submat((k-1),0,(k-1),(k-2)) = gamma.t();
  sigmamat_full((k-1),(k-1)) = 1;
  return sigmamat_full;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat subtract_y(arma::cube phi_arr,arma::vec ind, arma::mat y, int n){
  arma::mat subtracted_y(y.n_rows,y.n_cols);
  subtracted_y.row(0)=y.row(0);
  for (int j=1; j<n; j++){
    subtracted_y.row(j)=(y.row(j).t()-phi_arr.slice(ind(j))*y.row(j-1).t()).t();
  }
  return subtracted_y;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube subtract_b_cube(arma::cube b_cube, arma::cube phi_arr, arma::vec ind, int n){
  arma::cube subtracted_b_cube(size(b_cube), arma::fill::zeros);
  subtracted_b_cube.slice(0)=b_cube.slice(0);
  for (int j=1; j<n; j++){
    subtracted_b_cube.slice(j)=b_cube.slice(j)-phi_arr.slice(ind(j))*b_cube.slice(j-1);
  }
  return subtracted_b_cube;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List generatemu_posterior(arma::mat mu_prior_precision, arma::cube precision_full_arr,
                          arma::vec ind, arma::vec mu_prior_mu, arma::cube subtracted_b_cube, arma::mat subtracted_y, int n){
  arma::mat mu_posterior_precision = mu_prior_precision;
  arma::vec mu_posterior_mu = mu_prior_precision*mu_prior_mu;
  for (int nn = 0; nn <	n; nn++) {
    mu_posterior_precision += subtracted_b_cube.slice(nn).t()*
      precision_full_arr.slice(ind(nn))*subtracted_b_cube.slice(nn);
    mu_posterior_mu+=subtracted_b_cube.slice(nn).t()*precision_full_arr.slice(ind(nn))*
      (subtracted_y.row(nn).t());
  }

  arma::mat mu_posterior_sigmamat=inv_sympd(mu_posterior_precision);
  mu_posterior_mu=mu_posterior_sigmamat*mu_posterior_mu;
  return List::create(Rcpp::Named("mu_posterior_sigmamat") = mu_posterior_sigmamat,
                      Rcpp::Named("mu_posterior_mu") = mu_posterior_mu);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat center_y(arma::mat subtracted_y, arma::cube subtracted_b_cube, arma::vec mu, int n){
  arma::mat centered_y(n, subtracted_y.n_cols);
  for(int nn=0; nn<n; nn++){
    centered_y.row(nn)=(subtracted_y.row(nn).t()-subtracted_b_cube.slice(nn)*mu).t();
  }
  return centered_y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List generategamma_posterior(arma::vec ind, arma::cube sigmamat_conditional_arr,
                             arma::mat gamma_prior_precision, arma::vec gamma_prior_mu, arma::mat centered_y, int h, int k){
  arma::mat sigmamat_conditional_present = sigmamat_conditional_arr.slice(h);
  arma::mat precision_conditional_present = inv_sympd(sigmamat_conditional_present);
  arma::uvec seq_vech = find(ind==h);
  arma::mat centered_y_h = centered_y.rows(seq_vech);
  arma::mat gamma_posterior_precision = as_scalar(sum(pow(centered_y_h.col(k-1),2)))*precision_conditional_present + gamma_prior_precision;
  arma::mat gamma_posterior_sigmamat = inv_sympd(gamma_posterior_precision);
  arma::vec gamma_posterior_mu = gamma_posterior_sigmamat*(precision_conditional_present*(centered_y_h.cols(0,(k-2)).t()*centered_y_h.col(k-1))
                                                             +gamma_prior_precision*gamma_prior_mu);
  return List::create(Rcpp::Named("gamma_posterior_sigmamat") = gamma_posterior_sigmamat,
                      Rcpp::Named("gamma_posterior_mu") = gamma_posterior_mu);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat diff_centered_y(arma::mat centered_y, arma::mat gamma_mat, arma::vec ind, int k, int n){
  arma::mat diffed_centered_y(n,k-1);
  for(int nn=0; nn<n; nn++)
    diffed_centered_y.row(nn)=centered_y.submat(nn,0,nn,k-2)-centered_y(nn,k-1)*gamma_mat.col(ind(nn)).t();
  return diffed_centered_y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_w(arma::mat y, arma::vec mu, arma::vec ind, int H, int k, int n){
  arma::mat mu_mat(k,H);
  for(int h=0; h<H; h++){
    mu_mat.col(h) = mu.subvec(h*k,((h+1)*k-1));
  }
  arma::mat w(n,k);
  for(int nn=0; nn<n; nn++){
    w.row(nn) = y.row(nn) - mu_mat.col(ind(nn)).t();
  }
  return w;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List generatephi_posterior(arma::vec ind1, arma::mat w, arma::cube precision_full_arr, arma::cube phi_arr,
                           arma::mat phi_prior_precision, arma::vec phi_prior_mu, int k, int h, int kk){
  arma::uvec seq_vech = find(ind1==h);
  arma::mat w_sub_1=w.rows(seq_vech-1);
  arma::mat w_sub=w.rows(seq_vech);

  arma::vec seq_k = arma::linspace<arma::vec>(0,(k-1),k);
  arma::uvec seq_kk = find(seq_k!=kk);

  arma::mat w_sub_centered=w_sub-w_sub_1.cols(seq_kk)*phi_arr.slice(h).cols(seq_kk).t();

  arma::mat phi_posterior_sigmamat = inv_sympd(dot(w_sub_1.col(kk),w_sub_1.col(kk))*precision_full_arr.slice(h) + phi_prior_precision);
  arma::vec phi_posterior_mu = phi_posterior_sigmamat*((precision_full_arr.slice(h)*w_sub_centered.t())*w_sub_1.col(kk)
                                                         +phi_prior_precision*phi_prior_mu);

  return List::create(Rcpp::Named("phi_posterior_sigmamat") = phi_posterior_sigmamat,
                      Rcpp::Named("phi_posterior_mu") = phi_posterior_mu);
}

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double get_bovera1(arma::vec ind, arma::mat mu_present_mat, arma::cube precision_full_arr,
//                    arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn){
//   arma::vec mu_con = mu_present_mat.col(ind(nn-1)-1);
//   double a = as_scalar(sphericaldata.row(nn-1)*precision_full_arr.slice(ind(nn-1)-1)*sphericaldata.row(nn-1).t()+
//                        sphericaldata.row(nn-1)*phi_arr.slice(ind(nn)-1).t()*precision_full_arr.slice(ind(nn)-1)*phi_arr.slice(ind(nn)-1)*sphericaldata.row(nn-1).t());
//   double b = as_scalar(sphericaldata.row(nn-1)*precision_full_arr.slice(ind(nn-1)-1)*mu_con+
//                        (sphericaldata.row(nn-1)*phi_arr.slice(ind(nn)-1).t()*precision_full_arr.slice(ind(nn)-1)
//                           *(phi_arr.slice(ind(nn)-1)*mu_present_mat.col(ind(nn-1)-1)-mu_present_mat.col(ind(nn)-1)+augmenteddata_raw.row(nn).t())));
//   double bovera = b/a;
//   return bovera;
// }
//
//
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double get_boverann(arma::vec ind, arma::mat mu_present_mat, arma::cube precision_full_arr,
//                     arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn){
//   arma::vec mu_con = mu_present_mat.col(ind(nn-1)-1)+
//     phi_arr.slice(ind(nn-1)-1)*(augmenteddata_raw.row(nn-2).t()-mu_present_mat.col(ind(nn-2)-1));
//   double a = as_scalar(sphericaldata.row(nn-1)*precision_full_arr.slice(ind(nn-1)-1)*sphericaldata.row(nn-1).t()+
//                        sphericaldata.row(nn-1)*phi_arr.slice(ind(nn)-1).t()*precision_full_arr.slice(ind(nn)-1)*phi_arr.slice(ind(nn)-1)*sphericaldata.row(nn-1).t());
//   double b = as_scalar(sphericaldata.row(nn-1)*precision_full_arr.slice(ind(nn-1)-1)*mu_con+
//                        (sphericaldata.row(nn-1)*phi_arr.slice(ind(nn)-1).t()*precision_full_arr.slice(ind(nn)-1)
//                           *(phi_arr.slice(ind(nn)-1)*mu_present_mat.col(ind(nn-1)-1)-mu_present_mat.col(ind(nn)-1)+augmenteddata_raw.row(nn).t())));
//   double bovera = b/a;
//   return bovera;
// }
//
//
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// double get_boveran(arma::vec ind, arma::mat mu_present_mat, arma::cube precision_full_arr,
//                    arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn){
//   arma::vec mu_con = mu_present_mat.col(ind(nn-1)-1)+
//     phi_arr.slice(ind(nn-1)-1)*(augmenteddata_raw.row(nn-2).t()-mu_present_mat.col(ind(nn-2)-1));
//   double a = as_scalar(sphericaldata.row(nn-1)*precision_full_arr.slice(ind(nn-1)-1)*sphericaldata.row(nn-1).t());
//   double b = as_scalar(sphericaldata.row(nn-1)*precision_full_arr.slice(ind(nn-1)-1)*mu_con);
//   double bovera = b/a;
//   return bovera;
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double gen_r(arma::vec ind, arma::mat mu_present_mat, arma::cube precision_full_arr,
             arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn,
             int n, double r, int k){
  arma::vec mu_con;
  if(nn==0){
    mu_con = mu_present_mat.col(ind(nn));
  } else mu_con = mu_present_mat.col(ind(nn))+
    phi_arr.slice(ind(nn))*(augmenteddata_raw.row(nn-1).t()-mu_present_mat.col(ind(nn-1)));

  double a; double b;
  if(nn==n-1){
    a = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*sphericaldata.row(nn).t());
    b = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*mu_con);
  } else {
    a = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*sphericaldata.row(nn).t()+
      sphericaldata.row(nn)*phi_arr.slice(ind(nn+1)).t()*precision_full_arr.slice(ind(nn+1))*phi_arr.slice(ind(nn+1))*sphericaldata.row(nn).t());
    b = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*mu_con+
      (sphericaldata.row(nn)*phi_arr.slice(ind(nn+1)).t()*precision_full_arr.slice(ind(nn+1))
         *(phi_arr.slice(ind(nn+1))*mu_present_mat.col(ind(nn))-mu_present_mat.col(ind(nn+1))+augmenteddata_raw.row(nn+1).t())));

  }

  double bovera = b/a;
  double uniform1 = arma::randu();
  double uniform2 = arma::randu();
  double v = uniform1*exp(-0.5*a*pow((r-bovera),2.0));
  double root = sqrt(-2*log(v)/a);
  double upper = bovera + root;
  double lower = bovera + std::max(-bovera,-root);
  r = pow((pow(upper,k)-pow(lower,k))*uniform2+pow(lower,k),1/(double)k);
  return r;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double gen_r_group(arma::vec ind, arma::vec ind1, arma::mat mu_present_mat, arma::cube precision_full_arr,
             arma::cube phi_arr, arma::mat sphericaldata, arma::mat augmenteddata_raw, int nn,
             int n, double r, int k){
  arma::vec mu_con;
  if(nn==0){
    mu_con = mu_present_mat.col(ind(nn));
  } else mu_con = mu_present_mat.col(ind(nn))+
    phi_arr.slice(ind1(nn))*(augmenteddata_raw.row(nn-1).t()-mu_present_mat.col(ind(nn-1)));

  double a; double b;
  if(nn==n-1){
    a = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*sphericaldata.row(nn).t());
    b = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*mu_con);
  } else {
    a = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*sphericaldata.row(nn).t()+
      sphericaldata.row(nn)*phi_arr.slice(ind1(nn+1)).t()*precision_full_arr.slice(ind(nn+1))*phi_arr.slice(ind1(nn+1))*sphericaldata.row(nn).t());
    b = as_scalar(sphericaldata.row(nn)*precision_full_arr.slice(ind(nn))*mu_con+
      (sphericaldata.row(nn)*phi_arr.slice(ind1(nn+1)).t()*precision_full_arr.slice(ind(nn+1))
         *(phi_arr.slice(ind1(nn+1))*mu_present_mat.col(ind(nn))-mu_present_mat.col(ind(nn+1))+augmenteddata_raw.row(nn+1).t())));

  }

  double bovera = b/a;
  double uniform1 = arma::randu();
  double uniform2 = arma::randu();
  double v = uniform1*exp(-0.5*a*pow((r-bovera),2.0));
  double root = sqrt(-2*log(v)/a);
  double upper = bovera + root;
  double lower = bovera + std::max(-bovera,-root);
  r = pow((pow(upper,k)-pow(lower,k))*uniform2+pow(lower,k),1/(double)k);
  return r;
}





