#' Function to import DNA sequences data in the format required by \code{BayesANT}
#'
#' @param fasta.file File in .fas or .fasta format containing the DNA sequences.
#'                   The correct format for the annotation in the file should start with the ID of the sequence, followed by
#'                   a space and the word \code{Root;}. For example,
#'                   \code{>COLFF973-13 Root;Arthropoda;Insecta;Coleoptera}
#'                   is an annotation of 3 levels for the reference sequence with ID \code{COLFF973-13}
#' @param rank       Numeric argument indicating the number of taxonomic ranks used to construct the library.
#'                   Counting starts after the \code{Root;} in the annotation names.
#'                   Default is \code{rank = NULL}, which automatically selects the lowest level in the loaded taxonomy (the maximum length of the annotations).
#' @param rank_names Names for the taxonomic ranks, eg. "Class" or "Species". Default is \code{rank_names = NULL}, which automatically labels
#'                   the selected ranks as \code{"Level1"} up to level \code{"Levelx"}, where \code{"x"} is the value of the parameter \code{rank}
#'
#' @return An object of class \code{c("data.frame", "BayesANT.data")}
#' @export
#'
#' @examples
#'
#' rank <- 6
#' rank_names <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
#' data <- read.BayesANT.data("sequences.fasta", rank = rank, rank_names = rank_names)
#'
read.BayesANT.data <- function(fasta.file, rank = NULL, rank_names = NULL){

  ## Step 1 - Load the data in fasta format with DNA sequences.
  ## The automatic format in which they are loaded is a named list where DNA sequences are
  ## loaded as strings and are capitalized.
  data_fasta <- seqinr::read.fasta(file = fasta.file, forceDNAtolower = F, as.string = T)

  ## Step 2 - create the Taxonomic library
  annot <- lapply(data_fasta, function(x) attr(x, "Annot"))
  annot <- gsub(">", "", annot)

  # Extract the sequences IDs
  IDs <- names(data_fasta)

  # Extract the taxonomy
  taxa <- stringr::str_split(annot, pattern = ";", simplify = T)

  # Exclude the first row (which contains the Original Name and the Root)
  taxa <- taxa[, -1]

  # Restrict to the rank desired
  if(is.null(rank)){
    ## Automatically select the lowest rank, which is the maximum length for the annotations.
    lv <- ncol(taxa)
  } else {
    lv <- rank
    if(rank > ncol(taxa)){
      stop(cat("Value for rank = ", rank, "exceeds the length for the taxonomic annotations in the library. Specify a value for rank lower or equal to ", ncol(taxa)))
    }
  }
  taxa <- taxa[, 1:lv]

  # Specify the rank names
  if(is.null(rank_names)){
    rank_names <- paste0("Level", c(1:lv))
  }
  colnames(taxa) <- rank_names[1:lv]

  # Step 3 - Create the BayesANT library
  data <- data.frame(cbind(taxa, "DNA" = unlist(lapply(data_fasta, function(x) stringr::str_to_upper(x)))))
  class(data) <- c("data.frame","BayesANT.data")
  return(data)
}

