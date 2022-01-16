#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesANT
#sample_Taxonomy <- function(columns_seq_j, rho, M, Marginalprobs){
#
#  # Compute the probabilities
#  predicted_probs <- c(compute_probs(M=M, y=c(columns_seq_j), marginals = Marginalprobs$marginal_prob))^rho#
#
#  # Renormalize the probabilities
#  Marginalprobs$predicted_probs <- predicted_probs/sum(predicted_probs)

  # Obtain the final prediction
#  sample_one_taxon(Marginalprobs)

#}



#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesANT
#' @export
sample_Taxonomy <- function(x, s, rho, n_samples){
  # Extract the relevant information
  M <- x$M
  Priorprobs <- x$Priorprobs
  nucl <- x$nucl
  type_location <- x$type_location
  cols_to_drop <- NULL

  # Process the sequence
  if(type_location == "single"){
    s <- stringr::str_split(s, "")[[1]]
    s[!s %in% nucl] <- "-"

    if (!is.null(cols_to_drop)) {
      s <- s[-cols_to_drop]
    }
  } else if (type_location == "pairs"){
    s <- c(extract_2mers(seqDNA = s, nucl))
  }

  # Find which nucleotide to select
  sel <- build_CountMatrix(t(as.matrix(s)), nucl = nucl)

  # Compute the probabilities
  predicted_probs <- c(compute_probs(M=M, y=c(sel),priors = Priorprobs$prior_prob))^rho

  # Renormalize the probabilities
  Priorprobs$predicted_probs <- predicted_probs/sum(predicted_probs)

  # Sample the values now
  sampled_out <- matrix(NA, nrow = n_samples, ncol = ncol(Priorprobs) - 2)
  for(j in 1:n_samples){
    sampled_out[j,]<-sample_one_taxon(Priorprobs)
  }
  colnames(sampled_out) <- colnames(Priorprobs)[1:(ncol(Priorprobs) - 2)]
  return(sampled_out)
}



