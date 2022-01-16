
############### LOGPREDICTIVE
logPredictive = function(s, nodes, Nucl_counts, hyperpriors, nucl, cols_to_drop = NULL){
  p = length(s)
  # Extract the sequence
  sel = build_CountMatrix(t(as.matrix(s)), nucl = nucl)==1
  # Inizialize
  logProb =  rep(0, length(nodes))

  for(i in 1:length(nodes)){
    index_phy = which(names(Nucl_counts) == nodes[i])
    index_prior = which(names(hyperpriors) == nodes[i])

    if(grepl("_new", nodes[i])){
      # Use the prior predictive
      #if(is.null(cols_to_drop)){
      #  logProb[i] = sum(log(hyperpriors[[index_prior]][sel]))
      #} else {
      #  logProb[i] = sum(log(hyperpriors[[index_prior]][sel])[-cols_to_drop])
      #}
      if(is.null(cols_to_drop)){
        m1 = hyperpriors[[index_prior]]
      } else {
        m1= hyperpriors[[index_prior]][, -cols_to_drop]
      }
      logProb[i] = sum(log(m1[sel]/colSums(m1)))
    } else {
      # Use the posterior predictive
      index = which(names(Nucl_counts) == nodes[i])
      if(is.null(cols_to_drop)){
        m1 = (hyperpriors[[index_prior]] + Nucl_counts[[index_phy]])
      } else {
        m1 = (hyperpriors[[index_prior]][,-cols_to_drop] + Nucl_counts[[index_phy]][,-cols_to_drop])
      }
      logProb[i] = sum(log(m1[sel]/colSums(m1)))
    }

  }
  names(logProb) = nodes
  return(logProb)
}

####### PREDICT THE TREE

predict_Taxonomy <- function(s, Nucl_counts, terminal_nodes, Marginalprobs,
                             hyperpriors, cols_to_drop = NULL, nucl = c("-", "A", "C", "G", "T"),
                             return_wholetree = FALSE, type_location = "single") {
  final_pred <- data.frame()
  #best_clust <- rep(NA, 4)

  # Extract the sequence
  if(type_location == "single"){
    s <- stringr::str_split(s, "")[[1]]
    s[!s %in% nucl] <- "-"

    if (!is.null(cols_to_drop)) {
      s <- s[-cols_to_drop]
    }
  } else if (type_location == "pairs"){
    s <- c(extract_2mers(seqDNA = s))
  }


  # Compute prediction probabilities
  nuclprob <- logPredictive(s=s, nodes = terminal_nodes, Nucl_counts = Nucl_counts,
                            hyperpriors = hyperpriors, cols_to_drop = cols_to_drop, nucl = nucl)

  predProbs <- nuclprob + Marginalprobs$marginal_prob

  # Retrieve the probability across the whole phylogenetic tree
  Marginalprobs$predicted_probs <- exp(predProbs - max(predProbs)) / sum(exp(predProbs - max(predProbs)))

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










