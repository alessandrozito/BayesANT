### Pror Probablitites under the Pitman_Yor model
get_Priorprobs_PY = function(data, PYpars){

  data = data[, -ncol(data)]
  column_names = colnames(data)
  Priorprobs = dplyr::distinct(data.frame(data))

  # Level 1
  new = rep(paste0(column_names[1], "_new"), length(column_names))
  names(new) = column_names
  df_new = data.frame(t(new))

  # All the others
  for(l in 2:(length(column_names))){

    for(node in unique(Priorprobs[, l-1])){
      name_new = paste(node, column_names[l], "new", sep = "_")
      if(l == 2){
        new = c(Priorprobs[Priorprobs[, l-1] == node, 1:(l-1)][1],rep(name_new, length(column_names) - l + 1))
      } else {
        new = c(purrr::as_vector(Priorprobs[Priorprobs[, l-1] == node, 1:(l-1)][1,]), rep(name_new, length(column_names) - l + 1))
      }
      names(new) = column_names
      df_new = rbind(df_new, data.frame(t(new)))
    }

  }

  Priorprobs = as.data.frame(rbind(Priorprobs, df_new))

  ##### NOW, compute the prior probabilities

  ## Level 1
  # Extract the parameters
  alpha = PYpars[[1]]$param[1]
  sigma = PYpars[[1]]$param[2]
  clusters = PYpars[[1]]$frequencies

  # Probability of a new cluster
  prob_new = log(alpha + length(clusters)*sigma) - log(sum(clusters) + alpha)
  names(prob_new) = paste0(column_names[1], "_new")

  # Probability of an old cluster
  probPY = log(clusters - sigma) - log(sum(clusters) + alpha)
  df_prob = data.frame(names(c(prob_new, probPY)), c(prob_new, probPY))
  nm = c(column_names[1], paste0(column_names[1], "_prob"))
  colnames(df_prob) = nm
  # Merge to the Phylogenetic dataset
  Priorprobs = merge(Priorprobs, df_prob, by = column_names[1])


  #### Other levels
  for(p in 2:length(PYpars)){

    df_prob = data.frame()
    alpha = PYpars[[p]]$param[1]
    sigma = PYpars[[p]]$param[2]

    for(i in 1:length(PYpars[[p]]$frequencies)){
      upper = names(PYpars[[p]]$frequencies)[i]
      clusters = PYpars[[p]]$frequencies[[i]]

      # Probability of a new cluster
      prob_new = log(alpha + length(clusters)*sigma) - log(sum(clusters) + alpha)
      names(prob_new) = paste(upper, column_names[p], "new", sep = "_")

      # Probability of an old cluster
      probPY = log(clusters - sigma) - log(sum(clusters) + alpha)
      nm = c(column_names[p], paste0(column_names[p], "_prob"))
      df_prob_level = data.frame(names(c(prob_new, probPY)), c(prob_new, probPY))
      colnames(df_prob_level) = nm
      df_prob = rbind(df_prob, df_prob_level)

    }

    Priorprobs =  dplyr::left_join(Priorprobs, df_prob, by = column_names[p])
  }
  Priorprobs[is.na(Priorprobs)] = 0 # Substitute null values with zero, because prob is equal to one

  Priorprobs$prior_prob =  rowSums(Priorprobs[, grepl("prob",colnames( Priorprobs))])

  return(Priorprobs)
}

get_Priorprobs_noNewSpecies = function(ref_seq, PYpars){

  ref_seq = ref_seq[, -ncol(ref_seq)]
  column_names = colnames(ref_seq)
  Priorprobs = dplyr::distinct(ref_seq)

  ##### NOW, compute the prior probabilities

  ## Level 1
  # Extract the parameters
  #alpha = PYpars[[1]]$param[1]
  #sigma = PYpars[[1]]$param[2]
  alpha = 0
  sigma = 0
  clusters = PYpars[[1]]$frequencies

  # Probability of a new cluster
  prob_new = log(alpha + length(clusters)*sigma) - log(sum(clusters) + alpha)
  names(prob_new) = paste0(column_names[1], "_new")

  # Probability of an old cluster
  probPY = log(clusters - sigma) - log(sum(clusters) + alpha)
  df_prob = data.frame(names(c(prob_new, probPY)), c(prob_new, probPY))
  nm = c(column_names[1], paste0(column_names[1], "_prob"))
  colnames(df_prob) = nm
  # Merge to the Phylogenetic dataset
  Priorprobs = merge(Priorprobs, df_prob, by = column_names[1])


  #### Other levels
  for(p in 2:length(PYpars)){

    df_prob = data.frame()
    #alpha = PYpars[[p]]$param[1]
    #sigma = PYpars[[p]]$param[2]

    for(i in 1:length(PYpars[[p]]$frequencies)){
      upper = names(PYpars[[p]]$frequencies)[i]
      clusters = PYpars[[p]]$frequencies[[i]]

      # Probability of a new cluster
      prob_new = log(alpha + length(clusters)*sigma) - log(sum(clusters) + alpha)
      names(prob_new) = paste(upper, column_names[p], "new", sep = "_")

      # Probability of an old cluster
      probPY = log(clusters - sigma) - log(sum(clusters) + alpha)
      nm = c(column_names[p], paste0(column_names[p], "_prob"))
      df_prob_level = data.frame(names(c(prob_new, probPY)), c(prob_new, probPY))
      colnames(df_prob_level) = nm
      df_prob = rbind(df_prob, df_prob_level)

    }

    Priorprobs =  dplyr::left_join(Priorprobs, df_prob, by = column_names[p])
  }
  Priorprobs[is.na(Priorprobs)] = 0 # Substitute null values with zero, because prob is equal to one

  Priorprobs$prior_prob =  rowSums(Priorprobs[, grepl("prob",colnames( Priorprobs))])

  return(Priorprobs)
}


