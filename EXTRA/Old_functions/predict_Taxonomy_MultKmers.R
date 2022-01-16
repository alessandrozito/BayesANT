
log_MultKmers <- function(x_tilde, x){
  lgamma(sum(x)) - sum(lgamma(x)) + sum(lgamma(x + x_tilde)) - lgamma(sum(x + x_tilde))
}

logPredictive_MultKmers <- function(s, k,  nodes, Kmer_counts, hyperpriors, adjust_weights = FALSE, cols_to_drop = NULL){
  # Extract the Kmers relative to the sequence
  Kmers <- kmer::kcount(x= ape::as.DNAbin(s), k = k)
  n_kmers <- ncol(Kmers)

  if(adjust_weights == TRUE){
    weight_seq <- floor(length(s[s!="-"])/k)/(length(s[s!="-"]) - k + 1)
    Kmers <- Kmers*weight_seq
  } else {
    weight_seq <- 1
  }

  if(!is.null(cols_to_drop)){
    Kmers <- Kmers[, -cols_to_drop]
  }

  logProb =  rep(0, length(nodes))

  for(i in 1:length(nodes)){
    index_phy = which(names(Kmer_counts) == nodes[i])
    index_prior = which(names(hyperpriors) == nodes[i])

    if(grepl("_new", nodes[i])){
      # Use the prior predictive
      m1 = hyperpriors[[index_prior]]
    } else {
      # Use the posterior predictive
      m1 = hyperpriors[[index_prior]] + colSums(Kmer_counts[[index_phy]]$Kmers)
    }

    if(!is.null(cols_to_drop)){
      m1 = m1[, -cols_to_drop]
    }

    # Compute the probability
    logProb[i] = log_MultKmers(x_tilde = Kmers, x = m1)

  }
  names(logProb) = nodes
  return(logProb)
}


#### Make the prediction
predict_Taxonomy_MultKmers <- function(s, k, rho, Kmer_counts, terminal_nodes, Marginalprobs,
                                   hyperpriors,  adjust_weights = TRUE, cols_to_drop = NULL, nucl = c("-", "A", "C", "G", "T"),
                                   return_wholetree = FALSE) {
  final_pred <- data.frame()

  # Extract the sequence
  s <- str_split(s, "")[[1]]
  s[!s %in% nucl] <- "-"

  # Compute prediction probabilities
  predProbs <- logPredictive_MultKmers(s=s, k=k, nodes=terminal_nodes,
                                   Kmer_counts=Kmer_counts,
                                   hyperpriors=hyperpriors,
                                   adjust_weights = adjust_weights,
                                   cols_to_drop = cols_to_drop) + Marginalprobs$marginal_prob

  # Retrieve the probability across the whole phylogenetic tree
  Marginalprobs$predicted_probs <- exp(predProbs - max(predProbs)) / sum(exp(predProbs - max(predProbs)))

  # Recalibrate the probabilties according to rho
  Marginalprobs$predicted_probs <- Marginalprobs$predicted_probs^rho / sum(Marginalprobs$predicted_probs^rho)

  depth = (ncol(Marginalprobs) - 2)/2
  prediction <- rep(NA, depth)
  prediction_names <- rep(NA, depth)

  for(i in 1:length(prediction)){
    if(i == 1){
      df = Marginalprobs
    } else {
      df = Marginalprobs[Marginalprobs[, i-1] == best, ]
    }

    prob <- tapply(df$predicted_probs, df[,i], sum)
    best <- names(which.max(prob))
    prediction[i] <- prob[which.max(prob)]
    prediction_names[i] <- best
  }

  if(return_wholetree == FALSE){
    return(c(unname(prediction_names), unname(prediction)))
  } else {
    return(list("prediction" = c(unname(prediction_names), unname(prediction)), "whole_tree" = Marginalprobs))
  }

}













