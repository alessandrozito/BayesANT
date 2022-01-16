Process_missing_branches <- function(ref_seq, verbose = TRUE) {

  ## This function takes the  reference sequence dataset and processes the missing branches/taxa
  ## withing the taxonomy. In particular, the rule is the following
  ## 1) If no missing values, then do not do anything
  ## 2) If a node has both missing leaves and fully observed ones, remove all missing values
  ## 3) If a node is linked to missing leaves only, select one DNA sequence at random and
  ##    and name it as *_imputed. Eg. Lepidoptera_Order_imputed all down the last level

  L <- ncol(ref_seq) - 1 # Number of levels
  level_names <- colnames(ref_seq)[1:L]
  imputed_counter <- 0
  # Identify the missing labels
  missing_taxa <- complete.cases(ref_seq)
  missing_taxa_id <- which(!missing_taxa)

  if (all(missing_taxa)) {
    if (verbose == TRUE) {
      cat("No missing branches detected in the taxonomy \n")
    }
  } else {
    ref_seq_incomplete <- ref_seq[!missing_taxa, 1:L]

    for (j in 1:nrow(ref_seq_incomplete)) {
      # Find the lowest level that is not missing
      m <- max(which(!is.na(ref_seq_incomplete[j, ])))
      obs_taxa <- ref_seq_incomplete[j, m]
      # Check if at the level below all branches are missing
      na_vals <- is.na(ref_seq[ref_seq[, m] %in% obs_taxa, m + 1])

      if (all(na_vals) == TRUE) {
        #browser()
        #print(ref_seq[ref_seq[, m] %in% obs_taxa, 1:L])
        imputed_counter <- imputed_counter + 1
        # Sample one DNA sequence at random and plug it in the taxonomy with a missing name
        missing_taxa_id <- which(!missing_taxa & ref_seq[, m] %in% obs_taxa)
        id <- missing_taxa_id[sample(x=1:length(missing_taxa_id), size = 1)]
        # Create a fake name for the taxon at lower levels
        ref_seq[id, (m + 1):L] <- paste0(paste0(obs_taxa, "_", level_names[(m + 1):L]), "_imputed")
        missing_taxa[id] <- TRUE
        #print(ref_seq[ref_seq[, m] %in% obs_taxa, 1:L])
      }
    }
    if (verbose == TRUE) {
      cat(paste0("Number of imputed DNA sequences by default: ", imputed_counter, " \n"))
    }
  }
  return(ref_seq)
}

#a=Process_missing_branches(ref_seq)

#a[grepl("_imputed", a$Level_4),1:L]

#ref_seq[ref_seq$Level_2 %in% "Sp2-169", 3:L] = NA


#a[,1:L]



