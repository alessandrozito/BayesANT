#' Run the clustering algorithm for the DNA sequences in the test set
#'
#' @param x A \code{BayesANT} object
#' @param DNA Vector of DNA sequences to cluster
#' @param niter number of iteration to
#' @param rho Parameter to adjust the probabilities
#'
#' @return
#' @export
clusterSequences <- function(x, DNA, niter = 20, rho=1, cores = 1){

  # Number of levels in the taxonomy
  L = (ncol(x$Marginalprobs) - 1)/2

  #### Process the DNA sequences according to the model
  doParallel::registerDoParallel(cores)
  ### Quickly build the DNA sequence representation under our clustering model
  DNA_cols <- foreach(s = DNA, .combine = "rbind") %dopar% {
    c(build_CountMatrix(seqDNA_mat = t(as.matrix(AlignedSeq_to_vector(s=s,
                      type_location = x$type_location, nucl = x$nucl, cols_to_drop = NULL)))))
  }

  print("Finished processing the DNA sequences")
  #### Initialize the Taxonomy state
  taxonomy_state <- list("M" = x$M,
                         "Marginalprobs" = x$Marginalprobs,
                         "Nucl_counts" = x$Nucl_counts,
                         "hyperparameters" = x$hyperparameters,
                         "terminal_nodes" = x$terminal_nodes,
                         "new_taxa_detected" = 0,
                         "PY_pars" = x$PYpars)
  rownames(taxonomy_state$M) = x$terminal_nodes
  print("Initialize the Taxonomy state")

  ## Initialize the list to store posterior samples
  Clustering_Samples = vector(mode = "list", length = niter)
  current_clust = data.frame(matrix(NA, nrow = length(DNA), ncol = 2*L))
  colnames(current_clust) = colnames(x$Marginalprobs)[1:(2*L)]

  for(r in 1:niter){
    print("Start iterating")
    print(paste0("Iteration: ", r))
    for(j in 1:length(DNA)){
      # Create a matrix with the DNA counts
      y = matrix(DNA_cols[j, ], nrow = length(x$nucl))

      # Iterate across the sequences
      print(j)

      ## Remove the item from the state if we are boyond the first iteration
      if(r > 1){
        taxon_to_remove = c(unlist(current_clust[j,1:L]))
        taxonomy_state <- remove_seq_from_taxonomy_state(y=y,
                                                         taxon_to_remove = taxon_to_remove,
                                                         taxonomy_state = taxonomy_state)
      }

      ## Compute the taxon probabilities and sample it
      sampled_taxon = sample_Taxonomy(columns_seq_j = DNA_cols[j, ],
                                rho = 0.2,
                                M = taxonomy_state$M,
                                Marginalprobs = taxonomy_state$Marginalprobs)

      ## Update the taxonomy state
      taxonomy_state <- add_seq_to_taxonomy_state(taxonomy_state = taxonomy_state,
                                                  y = y,
                                                  taxon_assignment = sampled_taxon[1:L])
      taxon_assignment = c(taxonomy_state$taxon_new_name, sampled_taxon[(L+1):(2*L)])
      taxonomy_state <- taxonomy_state$taxonomy_state
      ## Insert the current sequence in the dataset
      current_clust[j,] = taxon_assignment

    }

    Clustering_Samples[[r]] = current_clust
  }
  return(Clustering_Samples)
}


