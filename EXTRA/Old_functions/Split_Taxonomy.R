Split_Taxonomy <- function(size, taxonomy, rank){
  n <- which(colnames(taxonomy) == rank)
  index.test <- c()
  nodes <- unique(taxonomy[,n])
  while(length(index.test) < size){

    # Sample one order at random
    sampled.nodes<- sample(nodes, 1)
    # Sample one observation with that order
    id <- sample(which(taxonomy[,n] == sampled.nodes), 1)

    if(!id %in% index.test){
      index.test <- c(index.test, id)
    }
  }
  sort(index.test)
}
