#### Make the prediction
predict_Taxonomy_not_alinged <- function(s, k, rho, ParameterMatrix, Priorprobs,
                                         adjust_Kmer_length, nucl, return_probs, n_top_taxa) {
  final_pred <- data.frame()

  # Extract the sequence
  s <- stringr::str_split(s, "")[[1]]
  s[!s %in% nucl] <- "-"

  # Extract the Kmers relative to the sequence
  Kmers <- kmer::kcount(x = ape::as.DNAbin(s), k = k)
  n_kmers <- ncol(Kmers)

  # Adjust the weights
  if (adjust_Kmer_length == TRUE) {
    weight_seq <- floor(length(s[s != "-"]) / k) / (length(s[s != "-"]) - k + 1)
    Kmers <- Kmers * weight_seq
  } else {
    weight_seq <- 1
  }

  # Compute prediction probabilities
  # Recalibrate the probabilties according to rho
  predProbs <- c(compute_probs_MultKmers(ParameterMatrix, Priorprobs$prior_prob, Kmers))^rho
  Priorprobs$predicted_probs <- predProbs / sum(predProbs)

  # Obtain the prediction
  out <- obtain_prediction(Priorprobs)

  if (return_probs == FALSE) {
    return(out)
  } else {
    depth <- (ncol(Priorprobs) - 2) / 2
    data_probs <- Priorprobs[, 1:depth]
    data_probs$leaf_prob <- predProbs
    data_probs <- data_probs[order(data_probs$leaf_prob)]
    return(list(
      "prediction" = out,
      "n_top_taxa" = utils::head(data_probs, n_top_taxa)
    ))
  }
}
