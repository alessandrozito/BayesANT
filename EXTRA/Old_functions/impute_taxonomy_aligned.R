AlignedSeq_to_vector <- function(s, type_location, nucl, cols_to_drop) {
  if (type_location == "single") {
    s <- stringr::str_split(s, "")[[1]]
    s[!s %in% nucl] <- "-"

    if (!is.null(cols_to_drop)) {
      s <- s[-cols_to_drop]
    }
  } else if (type_location == "pairs") {
    s <- c(extract_2mers(seqDNA = s))
  }
  s
}


impute_taxonomy_aligned <- function(j, ref_seq_missing, rho, M, Marginalprobs, nucl,
                                    type_location, adjust_weights, cols_to_drop, method = "prediction", ...) {
  L <- ncol(ref_seq_missing) - 1
  ## Function to impute the taxonomy via aligned prediction
  s <- ref_seq_missing[j, ncol(ref_seq_missing)]
  s <- AlignedSeq_to_vector(s, type_location, nucl, cols_to_drop) # Change the sequence to a vector
  sel <- build_CountMatrix(t(as.matrix(s)), nucl = nucl) # Extract the columns in the Matrix M

  # Look at where the sequence is missing
  tree_branch <- unlist(ref_seq_missing[j, -ncol(ref_seq_missing)])
  last_level <- max(which(!is.na(tree_branch))) # Level that is not missing from the taxonomy.

  # Restrict to just the child originating from that node.
  node <- ref_seq_missing[j, last_level]
  index_childred <- which(Marginalprobs[, last_level] == node)
  df <- Marginalprobs[index_childred, ]

  # Compute the probabilities at the leafs connected to the last true node
  predicted_probs <- c(compute_probs(
    M = M[index_childred, ], y = c(sel),
    marginals = Marginalprobs$marginal_prob[index_childred]))^rho
  df$predicted_probs <- predicted_probs / sum(predicted_probs) # Readjust them to account for rho

  # Allocate the sequence
  if (method == "prediction") {
    imputation_via_prediction(df, L, last_level, tree_branch)
  } else if (method == "cluster") {
    NULL
  }
}










# ref_seq_missing = model$ref_seq_missing
# M = model$M
# Marginalprobs = model$Marginalprobs
# type_location = x$type_location
# adjust_weight = x$adjust_weights
# cols_to_drop = x$cols_to_drop
# nucl = x$nucl
