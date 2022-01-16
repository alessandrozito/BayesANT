# Subtract from the Pitman-Yor frequencies
remove_seq_from_PYtree <-function(PY_pars, taxon_to_remove){
  L = length(PY_pars)
  is_removed = rep(0, L)

  ########### First level
  tab = PY_pars[[1]]$frequencies
  idx = which(names(tab) == taxon_to_remove[1])
  tab[idx] = tab[idx] - 1

  # Check if then node has been removed
  if(tab[idx] == 0){
      tab = tab[tab>0]
      is_removed[1] = 1
  }
  PY_pars[[1]]$frequencies = tab

  ########### All the other levels
  for(l in 2:L){
    idx_up = which(names(PY_pars[[l]]$frequencies) == taxon_to_remove[l-1])
    tab = PY_pars[[l]]$frequencies[[idx_up]]

    idx = which(names(tab) == taxon_to_remove[l])
    tab[idx] = tab[idx] - 1

    # Check if it is zero
    if(tab[idx] == 0){
      is_removed[l] = 1
      if(sum(tab)==0){
        # Remove entirely the cluster from the list
        PY_pars[[l]]$frequencies = PY_pars[[l]]$frequencies[-idx_up]
      } else {
        tab = tab[tab>0]
        PY_pars[[l]]$frequencies[[idx_up]] = tab
      }
    } else {
      PY_pars[[l]]$frequencies[[idx_up]] = tab
    }
  }
  return(list("PY_pars" = PY_pars, "is_removed"= is_removed))
}

########### Level 1
#PY_pars = taxonomy_state$PY_pars

#PY_pars_upd = remove_seq_from_PYtree(PY_pars, taxon_to_remove)


#PY_pars$Level_1$frequencies
#PY_pars_upd$Level_1$frequencies


#PY_pars$Level_2$frequencies$`Sp1-15`
#PY_pars_upd$Level_2$frequencies$`Sp1-15`


#PY_pars$Level_4$frequencies$NewTaxon_7
#PY_pars_upd$Level_4$frequencies$NewTaxon_8











