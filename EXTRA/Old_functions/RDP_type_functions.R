#### RDP-type functions
create_bag_of_words <- function(s, k, use_package){
  # k = kmer
  # s = sequence
  if(use_package){
    ## Use the kmer pakage. This is useful for cases where alignment gaps or
    ## special characters are available
    kmers <-kcount(as.DNAbin(str_split(data$DNAseq[1], "")[[1]]), k = 8)
    bag <- colnames(kmers)[kmers>0]
  } else {
    windows <- c(1:k)
    N <- length(s)
    kmers <- rep(NA,N-k+1)
    for(j in 1:length(kmers)){
      sel <- 1:k + (j-1)
      kmers[j] <- stringr::str_c(s[sel], collapse = "")
    }
    bag <- sort(unique(kmers))
  }
  return(bag)
}

get_kmers_in_library <- function(DNA, k, use_package = F, cores = 1){
  registerDoParallel(cores)
  # Extract the DNA list of sequences
  DNA_list <- stringr::str_split(DNA, pattern = "")
  # Find the bag of words for each sequence
  bags <- foreach(i=1:length(DNA_list)) %dopar% {
    create_bag_of_words(s = DNA_list[[i]], k = k,use_package = use_package)
  }
  return(bags)
}

#### Now, we are ready to create the kernel
get_Kmer_counts_RDP <- function(labels, kmer_bag){
  TaxaNames <- unique(labels)
  Nucl_counts <- vector(mode = "list", length(TaxaNames))
  for(i in 1:length(TaxaNames)){
    ## Extract the kmer sequence counts
    sel <- which(labels == TaxaNames[i])
    Nucl_counts[[i]] <- table(unlist(kmer_bag[sel]))
  }
  names(Nucl_counts) = TaxaNames
  return(Nucl_counts)
}


tune_hyperparameters_RDP <- function(tree, kmer_bag, labels, eps = 0.1, verbose = T){

  ## Create the hyperparameter probabilities in RDP
  hyperpriors = vector(mode = "list", length = nrow(tree))
  cl_col = ncol(tree)
  names(hyperpriors) = tree[, cl_col]

  ## Level of reference
  higher_cl = unique(tree[, cl_col - 1])
  new_phylum =paste0(colnames(tree)[1], "_new")

  niter = 0
  n_clust = length(higher_cl)
  verbose_step = max(round(n_clust/10), 1)

  for(cl in higher_cl){
    niter = niter + 1
    #print(cl)
    if(niter %% verbose_step == 0 & verbose == TRUE){
      cat(paste0("Number of Taxa processed: ", niter, " out of ", n_clust, "\n"))
    }
    if(grepl("_new$", cl)){

      if(cl == new_phylum){
        clusters = tree[, cl_col]
      } else {
        # Find the parent level
        parent_level <- which(tree[which(tree[,cl_col-1] == cl), ] == cl)[1] - 1
        parent <- tree[, parent_level][tree[, cl_col -1] == cl]
        clusters = tree[, cl_col][tree[, parent_level] == parent]
      }
      #select = (names(Nucl_counts) %in% clusters)
      indexes = which(names(hyperpriors) == cl)

    } else {
      # Cluster not new. Condition on the level above.
      clusters = tree[, cl_col][tree[, cl_col-1] == cl]
      indexes = which(names(hyperpriors) %in% clusters)
      #select = (names(Nucl_counts) %in% clusters)
    }

    sel <- which(labels %in% clusters)
    n_seqs <- length(sel)
    WordCounts <- table(unlist(kmer_bag[sel]))
    WordProbs <- (WordCounts + eps)/(n_seqs + 2*eps)
    WordProbs <- c(WordProbs, "other" = (eps)/(n_seqs + 2*eps))

    for(i in indexes){hyperpriors[[i]] = WordProbs}
  }
  return(hyperpriors)

}

get_RDP_WordProbs <- function(terminal_nodes, Nucl_counts, hyperparameters, labels, verbose = T){

  n_clust <- length(terminal_nodes)
  WordProbs <- vector(mode = "list", length = n_clust)
  names(WordProbs) <- terminal_nodes

  # Names of the Kmer counts.
  nm <- names(Nucl_counts)
  verbose_step = max(round(n_clust/10), 1)

  for(i in 1:n_clust){
    if(i %% verbose_step == 0 & verbose == TRUE){
      cat(paste0("Number of Leafs processed: ", i, " out of ", n_clust, "\n"))
    }

    cl <- terminal_nodes[i]
    if(grepl("_new$", cl)){
      ## New label use only the prior predictive distribution
      WordProbs[[i]] = log(hyperparameters[[i]])
    } else {
      ## Adjust based on the Kmer counts
      WordCounts <- Nucl_counts[which(nm==cl)][[1]]
      LeafSize <- sum(labels == cl)
      Priors <- hyperparameters[which(terminal_nodes==cl)][[1]]

      # Which kmers are observed both in the prior and in the counts?
      ID_obs_Kmers <- names(Priors) %in% names(WordCounts)
      ## Observed words
      PostProbs <- (WordCounts + Priors[ID_obs_Kmers])/(LeafSize + 1)
      ## Not Observed words
      PostProbs <- c(PostProbs, Priors[!ID_obs_Kmers]/(LeafSize + 1))
      WordProbs[[i]] <- log(PostProbs)
    }
  }
  return(WordProbs)
}

compute_PredictionProb_RPD <- function(bag, TaxaProbs){
  ## Observed words
  obs_words <- bag %in% names(TaxaProbs)
  p <- sum(TaxaProbs[names(TaxaProbs) %in% bag]) +
    sum(!obs_words) * TaxaProbs[names(TaxaProbs)=="other"]
  return(unname(p))
}

predict_Taxonomy_RDP <- function(s, rho, k, K, Priorprobs, return_wholetree=FALSE){

    ## Load the sequence and process it
    s <- stringr::str_split(s, "")[[1]]

    # Find the bag of words
    bag <- create_bag_of_words(s, k, use_package = F)

    # Compute the probabilities
    predicted_probs <- (unlist(lapply(K, function(x) compute_PredictionProb_RPD(bag, x)))
                        + Priorprobs$prior_prob)
    predicted_probs <- exp(predicted_probs - max(predicted_probs))
    predicted_probs <- (predicted_probs/sum(predicted_probs))^rho
    # Renormalize the probabilities
    Priorprobs$predicted_probs <-predicted_probs/sum(predicted_probs)

    # Obtain the final prediction
    out <- obtain_prediction(Priorprobs)

    if(return_wholetree == FALSE){
      return(out)
    } else {
      return(list("prediction" = out, "whole_tree" = Priorprobs))
    }

  }





