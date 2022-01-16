#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesANT
predict_Taxonomy_v3 <- function(s, kmers, M, K, Marginalprobs, nucl, type_location, return_wholetree=FALSE, cols_to_drop=NULL){

  ## Load the sequence and process it
  if(type_location == "single"){
    s <- stringr::str_split(s, "")[[1]]
    s[!s %in% nucl] <- "-"

    if (!is.null(cols_to_drop)) {
      s <- s[-cols_to_drop]
    }
  } else if (type_location == "pairs"){
    s <- c(extract_2mers(seqDNA = s))
  }

  # Fin which nucleotide to select
  sel <- build_CountMatrix(t(as.matrix(s)), nucl = nucl)

  # Compute the kmers
  new_kmers <- kmer::kcount(x= ape::as.DNAbin(s), k = kmers)

  # Compute the probabilities
  Marginalprobs$predicted_probs <- c(compute_probs_nucl_and_kmers(M=M, K=K, y=c(sel), marginals = Marginalprobs$marginal_prob, new_kmers = new_kmers))

  # Aggregate them upward
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
