add_seq_to_PYtree <- function(PY_pars, taxon_assignment, taxon_new_name=NULL){
  ### This function adds one sequence to the PY tree structure in the taxonomy
  ### Taxon new name applies only to the branch for which new species have been detected
  is_new = grepl("_new$", taxon_assignment) # Check if something is new within the taxonomy
  L = length(taxon_assignment)
  if(is.null(taxon_new_name)){
    taxon_new_name = taxon_assignment
  }
  for(l in 1:L){

    if(l==1){
      # First level
      counts = PY_pars[[l]]$frequencies
      if(!is_new[1]){
        count_id = which(names(counts) == taxon_assignment[l])
        counts[count_id] = counts[count_id] + 1
        PY_pars[[l]]$frequencies = counts
      } else {
        # You identified an entirely new branch in the taxonomy. Need to create the branch at
        # every level of the Pitman-Yor tree.After this, you can stop
        counts = c(counts, 1)
        names(counts)[length(counts)] = taxon_new_name[1]
        PY_pars[[l]]$frequencies = counts
        for(j in 2:L){
          freq = PY_pars[[j]]$frequencies
          freq$temp = 1
          names(freq$temp) = taxon_new_name[j]
          names(freq)[names(freq) == "temp"] = taxon_new_name[j-1]
          PY_pars[[j]]$frequencies = freq
        }
        break # All the branch is new, no need to check the remaining counts
      }

    } else {
      parent_id = which(names(PY_pars[[l]]$frequencies) == taxon_assignment[l-1])
      counts = PY_pars[[l]]$frequencies[[parent_id]]
      if(!is_new[l]){
        # The branch is observed up to here. Just add 1 to the counts
        count_id = which(names(counts) == taxon_assignment[l])
        counts[count_id] = counts[count_id] + 1
        PY_pars[[l]]$frequencies[[parent_id]] = counts
      } else {
        # We have identified a new node within an existing branch. Need to insert this into levels below
        counts = c(counts, 1)
        names(counts)[length(counts)] = taxon_new_name[l]
        PY_pars[[l]]$frequencies[[parent_id]] = counts
        if(l!=L){
          for(j in (l+1):L){
            freq = PY_pars[[j]]$frequencies
            freq$temp = 1
            names(freq$temp) = taxon_new_name[j]
            names(freq)[names(freq) == "temp"] = taxon_new_name[j-1]
            PY_pars[[j]]$frequencies = freq
          }
        }
        break  # All the branches below are new. Can stop here
      }
    }
  }
  return(PY_pars)
}

#PY_pars_new = add_seq_to_PYtree(PY_pars, taxon_assignment, taxon_new_name)




