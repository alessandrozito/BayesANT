load_FinBOL <- function(type_split = "sequence", level = "species"){
  if(type_split =="sequence"){
    name <- "onseqs_"
  } else if (type_split =="taxa"){
    name <- "ontaxa_"
  }

  "../TaxonomicClassifier/FinBOL/data/split_on_sequence/train_onseqs_species.fasta"
  data_train <- seqinr::read.fasta(paste0("../TaxonomicClassifier/FinBOL/data/split_on_",
                                  type_split,"/train_", name, level, ".fasta"))
  data_test <- seqinr::read.fasta(paste0("../TaxonomicClassifier/FinBOL/data/split_on_",
                                 type_split,"/test_", name, level, ".fasta"))

  ## Let's test our model on every scenario
  #### Process the data in the BayesANT format
  ref_seq <- data.frame(cbind(do.call("rbind", stringr::str_split(names(data_train), "[|]"))),
                        unlist(unname(lapply(data_train, function(x) stringr::str_to_upper(stringr::str_c(x, sep = "", collapse = ""))))))

  ref_seq_test <- data.frame(cbind(do.call("rbind", stringr::str_split(names(data_test), "[|]"))),
                             unlist(unname(lapply(data_test, function(x) stringr::str_to_upper(stringr::str_c(x, sep = "", collapse = ""))))))

  if(level == "species"){
    colnames(ref_seq) = colnames(ref_seq_test) = c("ID","Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "Species", "DNA")
  } else if(level == "genus"){
    colnames(ref_seq) = colnames(ref_seq_test) = c("ID","Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "DNA")
  }
  return(list("ref_seq" = ref_seq, "ref_seq_test"= ref_seq_test))
}
