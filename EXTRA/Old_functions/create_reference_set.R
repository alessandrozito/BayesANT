create_reference_set <- function(data_fasta, rank, rank_names){
  # Step 1 - Extract the taxonomic annotations.

  # Notice that the format of the file must be respected. See the guidelines
  annot <- lapply(data_fasta, function(x) attr(x, "Annot"))
  annot <- gsub(">", "", annot)

  # Extract the sequences IDs
  IDs <- names(data_fasta)

  # Extract the taxonomy
  taxa <- str_split(annot, pattern = ";", simplify = T)

  # Exclude the first row (which contains the Original Name and the Root)
  taxa <- taxa[, -1]

  # Restrict to the rank desired
  lv <- which(str_to_lower(rank) == str_to_lower(rank_names))
  if(length(lv)==0) {
    stop(cat("Value for rank = ", rank,"is not available among rank_names. Specify a value within rank names: c(", stringr::str_c(rank_names, collapse = ", "), ") or change rank_names"))
  }

  taxa <- taxa[, 1:lv]
  colnames(taxa) <- rank_names[1:lv]

  # Create the reference sequence set
  ref_seq <- data.frame(cbind("ID" = IDs, taxa, "DNA" = unlist(data_fasta)))

  return(ref_seq)
}
