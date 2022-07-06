#' Predict the taxonomic annotation for query DNA sequences via a BayesANT model.
#'
#' @param object an object of class \code{BayesANT}
#' @param rho Temperature parameter to re-calibrate the leaf probabilities.
#'            Default set to \code{rho = 0.1}.
#' @param return_probs Whether to return the first \code{n_top_taxa} leafs with
#'                     the highest probability probabilities. Default is
#'                     \code{return_probs = FALSE}.
#' @param n_top_taxa Number of leafs to return for the prediction. Default is
#'                   \code{n_top_taxa = 5}. Valid only if
#'                   \code{return_probs = TRUE}
#' @param cores Optional. Specify the number of cores for predicting DNA
#'              sequences in parallel. Relies on a combination of the packages
#'              \code{foreach} and \code{doParallel}. Default is
#'              \code{cores = 1}, which corresponds to a standard \code{for}
#'              loop.
#' @param verbose Monitor the number of sequences predicted.
#' @param ... Additional parameters
#' @param DNA DNA sequence to predict. Consider loading it with the function
#'            \code{read.BayesANT.testDNA}
#'
#' @return An object of class \code{data.frame} or \code{list}.
#' @importFrom foreach %dopar% foreach
#' @export
predict.BayesANT <- function(object,
                             DNA,
                             rho = 0.1,
                             return_probs = FALSE,
                             n_top_taxa = 5,
                             cores = 1,
                             verbose = TRUE,
                             ...) {

  # Error messages.
  stopifnot(
    class(object) == "BayesANT",
    is.numeric(rho),
    rho >= 0 & rho <= 1,
    is.logical(return_probs),
    is.logical(verbose),
    is.numeric(cores),
    is.numeric(n_top_taxa),
    is.character(DNA)
  )

  ## Register the number of cores in doParallel.
  ## Default is cores = 1, which is equivalent to a standard for loop
  doParallel::registerDoParallel(cores)

  ## Monitor the output
  seq_indexes <- c(1:length(DNA))
  tot <- length(DNA)
  if (tot >= 10) {
    verbose_step <- round(tot / 10)
  } else {
    verbose_step <- 9
  }


  # Begin prediction
  if (object$typeseq == "aligned") {
    # Verify that the length of the DNA is equal to the length of the
    # training sequences
    length_DNA <- unique(stringr::str_length(DNA))
    if (length(length_DNA) != 1) {
      stop("Query DNA sequences of different lengths are not allowed when
           typeseq = 'aligned' in the model.")
    }
    if (length_DNA != object$sequences_length) {
      stop(cat(
        "Query DNA sequences must be of length equal to ",
        object$sequences_length, "base pairs. \n"
      ))
    }

    i <- 1
    out <- foreach(i = seq_indexes) %dopar% {

      #  Print option valid only for cores = 1
      if (i %% verbose_step == 0 & verbose == TRUE) {
        cat("Number of sequences predicted = ", i, "/", tot, " \n")
      }
      predict_Taxonomy_aligned(DNA[i],
        rho = rho,
        ParameterMatrix = object$ParameterMatrix,
        Priorprobs = object$Priorprobs,
        type_location = object$type_location,
        nucl = object$nucl,
        return_probs = return_probs,
        n_top_taxa = n_top_taxa
      )
    }
  } else if (object$typeseq == "not aligned") {
    out <- foreach(i = seq_indexes) %dopar% {

      ##  Print option valid only for cores = 1
      if (i %% verbose_step == 0 & verbose == TRUE) {
        cat("Number of sequences predicted = ", i, "/", tot, " \n")
      }
      predict_Taxonomy_not_alinged(DNA[i],
        k = object$kmers,
        rho = rho,
        ParameterMatrix = object$ParameterMatrix,
        Priorprobs = object$Priorprobs,
        adjust_Kmer_length = object$adjust_Kmer_length,
        return_probs = return_probs,
        n_top_taxa = n_top_taxa,
        nucl = object$nucl
      )
    }
  }


  ##### Reformat the output
  if (return_probs == T) {
    # Create a separate data with the results, and a list containing the
    # n_top_taxa predictions
    df_results <- data.frame(do.call("rbind", lapply(out, function(x) {
      x$prediction
    })))
    top_n_probs <- lapply(out, function(x) x$n_top_taxa)
    names(top_n_probs) <- names(DNA)
  } else {
    df_results <- data.frame(do.call("rbind", out))
  }

  # Best prediction dataset
  levels <- object$level_names
  colnames(df_results) <- c(levels, paste0("Prob_", levels))
  prob_cols <- which(grepl("Prob_", colnames(df_results)))
  df_results[, prob_cols] <- sapply(df_results[, prob_cols], as.numeric)
  rownames(df_results) <- names(DNA)

  if (return_probs == T) {
    return(list("prediction" = df_results, "top_n_probs" = top_n_probs))
  } else {
    return(df_results)
  }
}
