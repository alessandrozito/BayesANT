impute_taxonomy_raw <- function(j, ref_seq_missing, rho, k, K, Marginalprobs,
                                nucl, adjust_weights, cols_to_drop, method ="prediction", ...){
  L = ncol(ref_seq_missing) - 1
  ## Function to impute the taxonomy using the kmer decomposition
  s <- ref_seq_missing[j, L+1]
  s <- stringr::str_split(s, "")[[1]]
  s[!s %in% nucl] <- "-"

  # Extract the Kmers relative to the sequence
  Kmers <- kmer::kcount(x= ape::as.DNAbin(s), k = k)
  n_kmers <- ncol(Kmers)

  # Adjust the weights
  if(adjust_weights == TRUE){
    weight_seq <- floor(length(s[s!="-"])/k)/(length(s[s!="-"]) - k + 1)
    Kmers <- Kmers*weight_seq
  } else {
    weight_seq <- 1
  }

  # Drop eventual columns
  if(!is.null(cols_to_drop)){
    Kmers <- Kmers[, -cols_to_drop]
  }

  # Look at where the sequence is missing
  tree_branch = unlist(ref_seq_missing[j, -(L+1)])
  last_level = max(which(!is.na(tree_branch))) # Level that is not missing from the taxonomy.

  # Restrict to just the child originating from that node.
  node = ref_seq_missing[j, last_level]
  index_childred = which(Marginalprobs[, last_level] == node)
  df = Marginalprobs[index_childred,]

  # Compute the probabilities at the leafs connected to the last true node
  predicted_probs <- c(compute_probs_MultKmers(K[index_childred, ],
                                               Marginalprobs$marginal_prob[index_childred],
                                               Kmers))^rho
  df$predicted_probs <- predicted_probs/sum(predicted_probs)  # Readjust them to account for rho

  # Allocate the sequence
  if(method=="prediction"){
    imputation_via_prediction(df, L, last_level, tree_branch)
  } else if(method=="cluster"){
    NULL
  }

}










#ref_seq_missing = model$ref_seq_missing
#M = model$M
#Marginalprobs = model$Marginalprobs
#type_location = x$type_location
#adjust_weight = x$adjust_weights
#cols_to_drop = x$cols_to_drop
#nucl = x$nucl


