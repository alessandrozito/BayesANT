build_CountMatrix = function(seqDNA_mat, nucl = c("-", "A", "C", "G", "T"), case_upper = TRUE){

  if(case_upper == FALSE){
    nucl = tolower(nucl)
  }
  CountMatrix = matrix(NA, nrow = length(nucl), ncol = ncol(seqDNA_mat))
  for(d in 1:length(nucl)){
    CountMatrix[d,] = colSums(seqDNA_mat == nucl[d])
  }
  return(CountMatrix)
}


