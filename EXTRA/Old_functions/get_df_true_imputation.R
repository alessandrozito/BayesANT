get_df_true_imputation <- function(ref_seq_full, ref_seq_processed){

  # Return the true value
  df_true = ref_seq_full
  id_to_save <- !complete.cases(ref_seq_processed)
  # Get the names of the levels
  level_names = names(ref_seq_full)[-ncol(ref_seq_full)]
  max_levels = length(level_names)

  for(l in 2:max_levels){
    # Find the indexes of imputed species
    ids = grepl("*_imputed", ref_seq_processed[,l])
    if(sum(ids)>0){
      for(j in which(ids)){

        true_taxa <- ref_seq_full[j,l]
        assigned_name <- ref_seq_processed[j,l]
        df_true[ref_seq_full[,l] == true_taxa, l] <- assigned_name
        higher_node = unique(ref_seq_full[ref_seq_full[,l]==true_taxa, l-1])

        all_clusters <- unique(ref_seq_full[ref_seq_full[,l-1]==higher_node, l])

        if(min(all_clusters%in%true_taxa)==0){
          other_taxa = all_clusters[all_clusters!=true_taxa]
          df_true[ref_seq_full[,l] %in% other_taxa, l] = paste0(higher_node, "_", level_names[l], "_new")
        }
      }
    }
  }
  return(df_true[id_to_save,1:max_levels])
}
