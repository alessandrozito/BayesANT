build_ParameterMatrix <- function(nodes, Nucl_counts, hyperparameters){
  M <- matrix(0, nrow = length(nodes), ncol = length(Nucl_counts[[1]]))

  for(i in 1:length(nodes)){
    index_phy = which(names(Nucl_counts) == nodes[i])
    index_prior = which(names(hyperparameters) == nodes[i])

    if(grepl("_new$", nodes[i])){
      m1 = hyperparameters[[index_prior]]
    } else {
      m1 <- hyperparameters[[index_prior]] + Nucl_counts[[index_phy]]
    }
    M[i,] <- c(t(t(m1)/colSums(m1)))
  }
  return(M)
}

build_ParameterMatrix_from_Kmers <- function(nodes, Kmer_counts, hyperparameters){
  K <-matrix(0, nrow = length(nodes), ncol = ncol(Kmer_counts[[1]]$Kmers))

  for(i in 1:length(nodes)){
    index_phy = which(names(Kmer_counts) == nodes[i])
    index_prior = which(names(hyperparameters) == nodes[i])

    if(grepl("_new$", nodes[i])){
      m1 = hyperparameters[[index_prior]]
    } else {
      m1 <- hyperparameters[[index_prior]] + colSums(Kmer_counts[[index_phy]]$Kmers)
    }
    K[i,] <- m1
  }
  return(K)
}

