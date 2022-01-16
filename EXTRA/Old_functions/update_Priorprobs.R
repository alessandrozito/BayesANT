
update_Priorprobs <- function(df_levels, PYpars){
  column_names = colnames(df_levels)
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
  # Merge to the levels dataset
  Priorprobs = dplyr::left_join(df_levels, df_prob, by = column_names[1])


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

  Priorprobs$prior_prob =  rowSums(Priorprobs[, grepl("prob",colnames(Priorprobs))])

  return(Priorprobs)
}




