load_Simulation_scenarios <- function(sc = 1){

  #### Load the files
  train_file <- paste0("../TaxonomicClassifier/Simulated_data/Ale_simulations/data/Scenario",
  sc, "/scenario", sc, "_train.fasta")
  test_file <- paste0("../TaxonomicClassifier/Simulated_data/Ale_simulations/data/Scenario",
                      sc, "/scenario", sc, "_test.fasta")
  train <- read.fasta(file = train_file)
  test <- read.fasta(file = test_file)

  #### Extract taxonomy
  taxa_train <- data.frame(stringr::str_split(names(train), "[|]", simplify = TRUE))
  taxa_test <- data.frame(stringr::str_split(names(test), "[|]", simplify = TRUE))
  colnames(taxa_train) <- colnames(taxa_test) <- c("ID", "Level1", "Level2", "Level3","Level4")

  ####  DNA
  taxa_train$DNA <- unname(unlist(lapply(train,function(x) stringr::str_c(stringr::str_to_upper(x),
                                                      collapse = ""))))
  taxa_test$DNA <- unname(unlist(lapply(test,function(x) stringr::str_c(stringr::str_to_upper(x),
                                                                          collapse = ""))))

  list("train" = taxa_train, "test" = taxa_test)

}
