#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec compute_probs(const arma::mat& M, const arma::vec& y,
                        const arma::vec& priors) {
  arma::vec logprobs = sum(log(M.cols( arma::find(y == 1) )), 1) + priors;
  double m = max(logprobs);
  arma::vec probs = exp(logprobs - m);
  return probs/sum(probs);
}

// [[Rcpp::export]]
double log_MultKmers_cpp(const arma::rowvec& x, const arma::rowvec& x_tilde){
  double logprob;
  logprob = lgamma(sum(x)) - sum(lgamma(x)) + sum(lgamma(x + x_tilde)) - lgamma(sum(x + x_tilde));
  return logprob;
}

// [[Rcpp::export]]
arma::vec compute_probs_MultKmers(const arma::mat& K, const arma::vec& priors,
                               const arma::rowvec& new_kmers){
  int n = K.n_rows;
  arma::vec logprobs_kmers(n);

  for(int i=0; i<n; i++){
    logprobs_kmers(i) = log_MultKmers_cpp(K.row(i), new_kmers) + priors(i);
  }
  double m = max(logprobs_kmers);

  // Compute the probabilities.
  logprobs_kmers = exp(logprobs_kmers - m);
  return logprobs_kmers/sum(logprobs_kmers);
}

//// [[Rcpp::export]]
//arma::vec compute_probs_nucl_and_kmers(const arma::mat& M, const arma::mat& K ,
//                                       const arma::vec& y,const arma::rowvec& new_kmers,
//                                       const arma::vec& priors){
  // probabilities from the multinomial nucleotide model
  //  arma::vec logprobs_nucl = sum(log(M.cols( arma::find(y == 1) )), 1);
  // probabilities from the kmer model
  //int n = logprobs_nucl.size();
  //arma::vec logprobs_kmers(size(logprobs_nucl));

  //for(int i=0; i<n; i++){
  //  logprobs_kmers(i) = log_MultKmers_cpp(K.row(i), new_kmers);
  //  }
  //arma::vec logprobs = logprobs_nucl + logprobs_kmers + priors;
  //double m = max(logprobs);
  //arma::vec probs = exp(logprobs - m);
  //return probs/sum(probs);
  //}





