get_df_true <- function(ref_seq, ref_seq_test){
  df_true = ref_seq_test
  
  # Get the names of the levels
  level_names = names(ref_seq)[-ncol(ref_seq)]
  max_levels = length(level_names)
  for(l in 1:max_levels){
    # Find the indexes of new species
    id = !(df_true[,l] %in% unique(ref_seq[,l]))
    n_replace = sum(id)
    if(n_replace>0){
      if(l == 1){
        name_clust = paste0(level_names[l], "_new")
        df_true[id, l:max_levels] = matrix(rep(name_clust, n_replace*max_levels), 
                                           nrow = n_replace, ncol = max_levels)
      } else {
        for(i in which(id==TRUE)){
          if(!grepl("_new", df_true[i, l-1])){
            name_clust = paste0(df_true[i, l-1], "_", level_names[l], "_new")
            df_true[i, l:max_levels] = rep(name_clust, max_levels-l + 1)
          }
        }
      }
    }
  }
  return(df_true)
}