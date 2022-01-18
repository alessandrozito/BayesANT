#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr %>%
#' @useDynLib BayesANT
predict_Taxonomy_aligned <- function(s, rho, ParameterMatrix, Priorprobs, nucl,
                                     type_location, return_probs, n_top_taxa) {

  ## Load the sequence and process it
  if (type_location == "single") {
    s <- stringr::str_split(s, "")[[1]]
    s[!s %in% nucl] <- "-"
  } else if (type_location == "pairs") {
    s <- c(extract_2mers(seqDNA = s, nucl))
  }

  # Find which nucleotide to select
  sel <- build_CountMatrix(t(as.matrix(s)), nucl = nucl)
  # Compute the probabilities
  # browser()
  predicted_probs <- c(compute_probs(M = ParameterMatrix, y = c(sel), priors = Priorprobs$prior_prob))^rho

  # Renormalize the probabilities
  Priorprobs$predicted_probs <- predicted_probs / sum(predicted_probs)

  # Obtain the final prediction
  out <- obtain_prediction(Priorprobs)

  if (return_probs == FALSE) {
    return(out)
  } else {
    depth <- (ncol(Priorprobs) - 2) / 2
    data_probs <- Priorprobs[, 1:depth]
    data_probs$leaf_prob <- predicted_probs
    data_probs <- data_probs[order(data_probs$leaf_prob)]

    return(list(
      "prediction" = out,
      "n_top_taxa" = utils::head(data_probs, n_top_taxa)
    ))
  }
}
