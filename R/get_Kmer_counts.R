create_gaps <- function(DNA, nucl) {
  DNA[!DNA %in% nucl] <- "-"
  return(DNA)
}

get_Kmer_counts <- function(k, nodes, DNAseq, nucl, adjust_Kmer_length = TRUE) {
  # Create the matrix of Kmer counts
  dnalist <- lapply(stringr::str_split(DNAseq, ""), function(x) create_gaps(x, nucl))

  Kmer_mat <- kmer::kcount(x = ape::as.DNAbin(dnalist), k = k)
  cols_zeros <- which(colSums(Kmer_mat) == 0)

  n_i <- unlist(lapply(dnalist, function(x) length(x[x != "-"])))

  if (adjust_Kmer_length == TRUE) {
    # Adjust for the likelihood contribution of each kmer

    weights <- floor(n_i / k) / (n_i - k + 1)
    Kmer_mat <- Kmer_mat * weights
  } else {
    weights <- rep(1, length(n_i))
  }

  # Extract the species names
  species_names <- sort(unique(nodes))
  Kmer_counts <- vector(mode = "list", length = length(species_names))
  for (i in 1:length(species_names)) {
    idx <- (nodes == species_names[i])
    if (sum(idx) == 1) {
      Kmer_counts[[i]]$Kmers <- t(as.matrix(Kmer_mat[idx, ]))
    } else {
      Kmer_counts[[i]]$Kmers <- Kmer_mat[idx, ]
    }
    Kmer_counts[[i]]$weights <- weights[idx]
  }
  names(Kmer_counts) <- species_names

  return(Kmer_counts)
}
