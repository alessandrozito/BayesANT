extract_2mers <- function(seqDNA, nucl, string = TRUE){
  if(string == TRUE){
    # Feeding a string format. Need to turn it into a vector
    seqDNA <- stringr::str_split(seqDNA, "")[[1]]
    seqDNA[!seqDNA %in% c("-", "A", "C", "G", "T")] <- "-"
  }
  seqlen <- length(seqDNA)
  seq2mer<-c(apply(cbind(seqDNA[1:(seqlen-1)], seqDNA[2:seqlen]), 1, function(x) paste0(x, collapse ="")))
  return(seq2mer)
}

# Vectorize it on the DNA sequence
extract_2mers <- Vectorize(extract_2mers, vectorize.args = "seqDNA", USE.NAMES = FALSE)

get_2mer_counts <- function(nodes, DNAseq, nucl){

  # Extract the unique families
  unique_nodes = sort(unique(nodes))

  # Initialize an empty list
  Nucl_counts = vector(mode = "list", length = length(unique_nodes))
  names(Nucl_counts) = unique_nodes

  for(i in 1:length(unique_nodes)){

    indexes = which(unique_nodes[i] == nodes)
    # Extract the 2mer matrix
    matDNA <- t(extract_2mers(DNAseq[indexes], nucl = nucl))
    # Build the count matrix
    Nucl_counts[[i]] = build_CountMatrix(matDNA, nucl = nucl)
  }
  return(Nucl_counts)
}




