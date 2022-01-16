add_seq_to_taxonomy_state <- function(y, taxon_assignment, taxonomy_state){

  L <- length(taxon_assignment)
  level_names = names(taxonomy_state$Marginalprobs)[1:L]
  is_new = grepl("_new$", taxon_assignment[1:L])

  taxon_new_name <- taxon_assignment
  ############################################################
  ### STEP 1 - Update the Matrix M of DNA sequence score
  ############################################################
  if(all(!is_new)){
    ## Easy job. Just need to increment by one the row of M corresponding to that matrix
    index_phy = which(names(taxonomy_state$Nucl_counts) == taxon_assignment[L])
    index_prior = which(names(taxonomy_state$hyperparameters) == taxon_assignment[L])
    # Update the Nucleotide counts
    taxonomy_state$Nucl_counts[[index_phy]] <- taxonomy_state$Nucl_counts[[index_phy]] + y
    # Update the matrix M
    index_M = which(rownames(taxonomy_state$M) == taxon_assignment[L])
    m1 <- taxonomy_state$hyperparameters[[index_prior]] + taxonomy_state$Nucl_counts[[index_phy]]
    taxonomy_state$M[index_M, ] <- c(t(t(m1)/colSums(m1)))

  } else {
    ## A new taxon has been detected. We need to add a new hyperprior,
    ## a new element in the nucleotide counts and a new and updated row in the matrix M

    ## Create a new Nucleotide counts kernel
    #id = L - sum(is_new)
    #if(id==0){
    #  name = level_names[1]
    #} else {
    #  name =  paste(level_names[id],taxon_assignment[id], sep = ':')
    #}
    # Name of the new cluster
    clust = paste0("NewTaxon_", taxonomy_state$new_taxa_detected)
    taxon_new_name[is_new] = clust
    # Add 1 to the state, to distinguish it from future clusters
    taxonomy_state$new_taxa_detected = taxonomy_state$new_taxa_detected+1
    # Add the new matrix to the Nucl count list
    taxonomy_state$Nucl_counts$temp = y
    names(taxonomy_state$Nucl_counts)[names(taxonomy_state$Nucl_counts)=="temp"] = clust

    ## Create now the new hyperparameters. Note: one for each level below the new node
    hyper_new = taxonomy_state$hyperparameters[[which(names(taxonomy_state$hyperparameters) == taxon_assignment[L])]]

    for(l in which(is_new)){
      if(l==L){
        # If the node is at the last level
        taxonomy_state$hyperparameters$temp = hyper_new
        names(taxonomy_state$hyperparameters)[names(taxonomy_state$hyperparameters) == "temp"] = clust
        # Append the new observation to the matrix M
        m1 <- hyper_new + y
        taxonomy_state$M <- rbind(taxonomy_state$M , c(t(t(m1)/colSums(m1))))
        rownames(taxonomy_state$M)[nrow(taxonomy_state$M)] <-clust
      } else {
        # If the node is not at the last level. This creates subnodes
        taxonomy_state$hyperparameters$temp = hyper_new
        new_clust_name = paste0(clust, "_", level_names[l+1], "_new")
        names(taxonomy_state$hyperparameters)[names(taxonomy_state$hyperparameters) == "temp"] = new_clust_name
        # Append the new observation to the matrix M
        taxonomy_state$M <- rbind(taxonomy_state$M , c(t(t(hyper_new)/colSums(hyper_new))))
        rownames(taxonomy_state$M)[nrow(taxonomy_state$M)] <- new_clust_name
      }
    }
  }

  ############################################################
  ### STEP 2 - Update the Marginal Probabilities
  ############################################################

  ## First of all, we need to add the new sequence to the PY counts tree
  taxonomy_state$PY_pars <- add_seq_to_PYtree(taxonomy_state$PY_pars, taxon_assignment, taxon_new_name)
  ## Then, we need to verify if we included new taxa
  df_levels = taxonomy_state$Marginalprobs[, 1:L]
  if(sum(is_new)>0){
    for(l in which(is_new)){
      if(l==L){
        # Nothing happens. Just need to add this vector to the marginal probs
        df_levels = rbind(df_levels, taxon_new_name)
      } else {
        # Include also the new potential clusters if the node is discovered at a higher level
        potential_new_node = taxon_new_name
        potential_new_node[(l+1):L] = paste0(clust, "_", level_names[l+1], "_new")
        df_levels = rbind(df_levels, potential_new_node)
      }
    }
  }
  taxonomy_state$Marginalprobs = update_Marginalprobs(df_levels, PYpars = taxonomy_state$PY_pars)

  ## Return the updated taxonomy state
  return(list("taxonomy_state" = taxonomy_state, "taxon_new_name" = taxon_new_name))
}





#a = update_Marginalprobs(df_levels, PYpars = taxonomy_state$PY_pars)
#b = get_Marginalprobs(ref_seq = rbind(ref_seq, c(taxon_new_name, 'a')), PYpars = taxonomy_state$PY_pars)
#all(a %>% dplyr::arrange(Level_4) ==b %>% dplyr::arrange(Level_4))




