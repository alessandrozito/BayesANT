
#### Function to mutate a sequence at random
mutate_sequence <- function(sequence, eps, nucleotides){
  mutations <- as.logical(matrix(rbinom(n = ncol(sequence) * nrow(sequence),
                                        prob = eps, size = 1),
                                 nrow = nrow(sequence), ncol = ncol(sequence)))
  nucleotides_to_mutate <- sequence[mutations]
  if (length(nucleotides_to_mutate) > 0) {
    for (i in 1:length(nucleotides_to_mutate)) {
      nc <- nucleotides_to_mutate[i]
      nucleotides_to_mutate[i] <- sample(x = nucleotides[nucleotides != nc],
                                         size = 1, replace = TRUE)
    }
    sequence[mutations] <- nucleotides_to_mutate
  }

  return(sequence)
}


#' Function to simulate a reference sequence library from a coalescent process
#'
#' @param size Total number of sequences to generate
#' @param alphas alphas parameters. The size must be equal to n_levels
#' @param sigmas sigmas parameters. The size must be equal to n_levels
#' @param n_levels Number of levels to construct
#' @param length_conserved_region Length of the conserved region
#' @param nucleotides Nucleotides for the simulation. Default is c('a', 'c', 'g', 't')
#' @param level_names Names of the levels
#' @param length_seq Length of the sequence
#' @param rates Rates of mutation in the coalescent process
#' @param rate_lastnode Rates of mutation in the coalescent process
#' @param verbose Print the state of the simulation
#' @param root Root sequence from which to start he simulation
#' @param seed Set the internal seed for the simulation
#'
#' @return
#' @examples
#'
#' size = 200
#' n_levels = 3
#' alphas = c(2,1, 1)
#' sigmas = c(0,0, 0)
#' Q = get_K80mat(0.5,0.5)
#' Q_list = list(Q, Q, Q)
#' length_conserved_region = 50
#' length_seq = 300
#' rate_lastnode = 0.03
#' ref_seq <- simulate_library(size, n_levels, alphas, sigmas,
#'                  length_seq, rates, rate_lastnode, length_conserved_region)
#'
#' @export
#'
simulate_library <- function(size, n_levels, alphas, sigmas, length_seq,
                                           rates, rate_lastnode, length_conserved_region = NULL,
                                           nucleotides = c("a", "c", "g", "t"), level_names = NULL,
                                           verbose = FALSE, root = NULL, seed = NULL) {


  ## Set the seed
  if(!is.null(seed)){
    set.seed(seed)
  }

  # STEP 0 - simulate the ancestral sequence
  if(is.null(root)){
    root <- sample(x = nucleotides, size = length_seq, replace = TRUE)
  }

  # Exclude the conserved region
  if (!is.null(length_conserved_region)) {
    conserved_region <- t(as.matrix(root[c(1:length_conserved_region)]))
    root <- root[-c(1:length_conserved_region)]
  }

  # Name the levels
  if (is.null(level_names)) {
    level_names <- paste0("Level_", c(1:n_levels))
  }

  # Initialize the reference sequence dataset
  ref_seq <- data.frame(matrix(rep(NA, size * n_levels), nrow = size, ncol = n_levels))
  colnames(ref_seq) <- level_names

  ####################################### Generate the core sequences

  for(l in 1:n_levels){
    if(verbose==TRUE){cat(paste0("Simulating level: ", l, "\n"))}

    # Step 1 - Simualte from the Pitman-Yor for the given level
    if(l == 1){
      # Do the first level
      counts <- samplePY(size = size, alpha = alphas[l], sigma = sigmas[l])
    } else {
      # Sample all the other levels
      counts = sapply(counts, function(x) samplePY(size = x, alpha = alphas[l], sigma = sigmas[l]))
    }
    # Branches to split the ref sequence

    if(l==1){
      n_branches = length(counts)
    } else {
      n_branches = lapply(counts, length)
    }

    # Step 2 - Simulate a coalescent process for every cluster
    if(l == 1){
      coal <- ape::rcoal(n_branches)
      sequence <- as.character(phangorn::simSeq(coal, type="DNA",
                                                rootseq = root,
                                                l = length(root),
                                                rate = rates[l]))
    } else {
      whole_seq <- NULL
      for(j in 1:length(n_branches)){
        if(n_branches[j] == 1){
          temp_mat <- t(as.matrix(sequence[j, ]))
        } else {
          coal <- ape::rcoal(n_branches[j])
          temp_mat <- as.character(phangorn::simSeq(coal, type="DNA",
                                                    rootseq = sequence[j, ],
                                                    l = length(root),
                                                    rate = rates[l]))
        }
        whole_seq <- rbind(whole_seq, temp_mat)
      }
      sequence <- whole_seq
    }


    # Unlist the counts and add to the node level
    counts <- unlist(counts)
    node_names <- rep(paste0("V*-", c(1:nrow(sequence)), ",", l), counts)
    ref_seq[, l] <- node_names

  }

  ####################################### Expand sequence to generate the full one
  #sequence <- sequence[rep(1:nrow(sequence), counts), ]
  #sequence <- mutate_sequence(sequence, eps = epsilon_lastnode, nucleotides = nucleotides)
  final_seq <- NULL
  for(j in 1:nrow(sequence)){
    if(counts[j] == 1){
      final_seq <- rbind(final_seq, sequence[j, ])
    } else {
      # Generate a final distribution, possibly with a low rate.
      coal <- ape::rcoal(n = counts[j])
      newmat <- as.character(phangorn::simSeq(coal, type="DNA",
                                              rootseq = sequence[j, ],
                                              l = length(root),
                                              rate = rate_lastnode))
      newmat <- mutate_sequence(newmat, eps = rate_lastnode, nucleotides = nucleotides)
      final_seq <- rbind(final_seq, newmat)
    }
  }

  # RESTORE THE CONSERVED REGION
  if (!is.null(length_conserved_region)) {
    final_seq <- cbind(conserved_region[rep(1, nrow(final_seq)), ], final_seq)
  }

  # Merge the nucleotide strings
  all_seq <- apply(final_seq, 1, function(x) stringr::str_c(x, collapse = ""))
  ref_seq$DNAseq <- stringr::str_to_upper(all_seq)
  return(ref_seq)

}




