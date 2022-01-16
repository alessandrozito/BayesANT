remove_seq_from_taxonomy_state <- function(y, taxon_to_remove, taxonomy_state){
  L = length(taxon_to_remove)

  #### Step 1 - remove the item from the Pitman-Yor tree
  PY_pars_upd <- remove_seq_from_PYtree(taxonomy_state$PY_pars, taxon_to_remove)
  taxonomy_state$PY_pars <- PY_pars_upd$PY_pars

  #### Step 2 - Downgrade the Marginalprobs
  df_levels = taxonomy_state$Marginalprobs[,1:L]
  if(any(PY_pars_upd$is_removed==1)){
    lines_to_remove = grep(paste(taxon_to_remove[PY_pars_upd$is_removed==1], collapse = "|"), df_levels[,L])
    if(length(lines_to_remove)!= sum(PY_pars_upd$is_removed)){
      print("Error in nodes removed")
      break
    }
    df_levels = df_levels[-lines_to_remove,]
  }
  taxonomy_state$Marginalprobs = update_Marginalprobs(df_levels, taxonomy_state$PY_pars)


  #### Step 3 - Remove the Nucleotide counts
  index_nucl = which(names(taxonomy_state$Nucl_counts) == taxon_to_remove[L])
  if(PY_pars_upd$is_removed[L] == 0){
    # Subtract by one
      taxonomy_state$Nucl_counts[[index_nucl]] = taxonomy_state$Nucl_counts[[index_nucl]] - y
  } else {
    # Remove entirely
    taxonomy_state$Nucl_counts = taxonomy_state$Nucl_counts[-index_nucl]
  }

  #### Step 4 - Remove the hyperparameter (only if the cluster is removed)
  if(any(PY_pars_upd$is_removed==1)){
    idx = which(! names(taxonomy_state$hyperparameters) %in% taxonomy_state$Marginalprobs[,L])
    taxonomy_state$hyperparameters = taxonomy_state$hyperparameters[-idx]
  }

  #### Step 5 - Downgrade the Matrix M
  if(any(PY_pars_upd$is_removed==1)){
    idx = which(! rownames(taxonomy_state$M) %in% taxonomy_state$Marginalprobs[,L])
    taxonomy_state$M = taxonomy_state$M[-idx,]
  } else {
    index_phy = which(names(taxonomy_state$Nucl_counts) == taxon_to_remove[L])
    index_prior = which(names(taxonomy_state$hyperparameters) == taxon_to_remove[L])
    # Update the matrix M
    index_M = which(rownames(taxonomy_state$M) == taxon_to_remove[L])
    m1 <- taxonomy_state$hyperparameters[[index_prior]] + taxonomy_state$Nucl_counts[[index_phy]]
    taxonomy_state$M[index_M, ] <- c(t(t(m1)/colSums(m1)))
  }
  return(taxonomy_state)
}
