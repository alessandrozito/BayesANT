#' Load query DNA sequences to be used by the function \code{predict.BayesANT}
#'
#' @param fasta.file File in \code{.fasta} or \code{.fas}
#'                   format containing the DNA sequences.
#' @param ... Additional parameters.
#'
#' @return A named vector of strings
#' @export
read.BayesANT.testDNA <- function(fasta.file, ...) {
  # Load the test functions to be predicted by BayesANT
  data_fasta <- seqinr::read.fasta(
    file = fasta.file, forceDNAtolower = F,
    as.string = T
  )
  DNA <- unlist(lapply(data_fasta, function(x) stringr::str_to_upper(x)))
  names(DNA) <- names(data_fasta)
  return(DNA)
}
