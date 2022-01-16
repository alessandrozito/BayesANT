load_its2 <- function(split, level){
  # Load its2
  its2 <- readRDS("../TaxonomicClassifier/its2/data/its2.rds.gzip")
  file_train <- paste0("../TaxonomicClassifier/its2/data/its2-", split,"-train-test/its2-",level,"-train.txt")
  file_test <- paste0("../TaxonomicClassifier/its2/data/its2-", split,"-train-test/its2-",level,"-test.txt")


  trainIDs <- suppressMessages(readr::read_table(file_train, col_names = "ID"))
  testIDs <- suppressMessages(readr::read_tsv(file_test, col_names = "ID"))

  df_train <-merge(trainIDs, its2, by = "ID", sort = FALSE)
  df_test <-merge(testIDs, its2, by = "ID", sort = FALSE)

  # Ref seq
  ref_seq <- df_train[, -c(1,2,3)]
  ref_seq_test <- df_test[, -c(1,2,3)]
  col <- which(colnames(ref_seq)==stringr::str_to_title(level))
  ref_seq <- ref_seq[, c(1:col, ncol(ref_seq))]
  ref_seq_test <- ref_seq_test[, c(1:col, ncol(ref_seq_test))]


  if(length(colnames(ref_seq)) > 4){
    # TRAIN
    n_orders <- rep(1, length(unique(ref_seq$Family)))
    fams <- unique(ref_seq$Family)
    for(i in 1:length(fams)){
      n_orders[i] <- length(unique((ref_seq %>% filter(Family==fams[i]))$Order))
    }

    names(n_orders) = fams
    names_double_family <- names(which(n_orders > 1))
    for(j in 1:length(names_double_family)){
      ids <- which(ref_seq$Family == names_double_family[j])
      ref_seq$Family[ids] <-  paste0(ref_seq$Family[ids], rep("_from_",length(ids)),  ref_seq$Order[ids])
    }
    # TEST
    n_orders <- rep(1, length(unique(ref_seq_test$Family)))
    fams <- unique(ref_seq_test$Family)
    for(i in 1:length(fams)){
      n_orders[i] <- length(unique((ref_seq_test %>% filter(Family==fams[i]))$Order))
    }

    names(n_orders) = fams
    names_double_family <- names(which(n_orders > 1))
    for(j in 1:length(names_double_family)){
      ids <- which(ref_seq_test$Family == names_double_family[j])
      ref_seq_test$Family[ids] <-  paste0(ref_seq_test$Family[ids], rep("_from_",length(ids)),  ref_seq_test$Order[ids])
    }
  }

  ## One sequence is really problematic.. need to remove it (only has NNNNNNNACG...)
  if("SH016704.07FU_JQ761044" %in% trainIDs$ID){
    rm_id <- which(trainIDs$ID == "SH016704.07FU_JQ761044")
    trainIDs <-trainIDs[-rm_id,]
    ref_seq <- ref_seq[-rm_id, ]
  } else if ("SH016704.07FU_JQ761044" %in% testIDs$ID){
    rm_id <- which(testIDs$ID == "SH016704.07FU_JQ761044")
    testIDs <-testIDs[-rm_id,]
    ref_seq_test <- ref_seq_test[-rm_id, ]
  }
  return(list("ref_seq" = ref_seq, "ref_seq_test" = ref_seq_test, "train_id" = trainIDs$ID, "test_id" = testIDs$ID))
}
