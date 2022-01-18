#' Summary function for a model of class \code{BayesANT}
#'
#' @param object An object of class \code{BayesANT}
#' @param ... Additional parameters
#'
#' @details A summary of the model details, including parameter estimation
#' @export
summary.BayesANT <- function(object, ...){
  ## Create model summary

  # Create summary table for the taxonomy and the estimates
  res <- NULL
  tot_nodes <- 0
  tot_new <- 0
  for(j in 1:(ncol(object$data) - 1)){
    nodes <- unique(object$Priorprobs[,j])
    n_observed <- sum(!grepl("_new$", nodes))
    n_new <- length(nodes) - n_observed
    #Update new nodes
    tot_nodes <-tot_nodes + length(nodes)
    tot_new <-tot_new + n_new

    if(object$newtaxa == TRUE){
      alpha <- round(object$PYpars[[j]]$param[1], 2)
      sigma <- round(object$PYpars[[j]]$param[2], 2)
      loglik <- round(object$PYpars[[j]]$loglik, 2)
    } else {
      alpha <- sigma <- loglik <- "."
    }

    res <-rbind(res, c(n_observed, n_new, alpha, sigma, loglik))
  }
  colnames(res) <- c("Observed taxa", "New Taxa", "alpha", "sigma", "loglik")
  rownames(res) <- object$level_names

  # Look at the DNA sequences
  if(object$typeseq == "aligned"){
    lDNA <- stringr::str_length(object$data$DNA[1])
  } else {
    lDNA <- mean(stringr::str_length(object$data$DNA))
  }

  # Print now the final message
  cat(paste("BayesANT classifier", paste(rep(xml2::xml_text(xml2::read_html(paste0("<x>", "&#128028;", "</x>"))), 4), collapse = "")),
      paste0("\nNumber of DNA sequences in the library: ", nrow(object$data)),
      paste0("Type of sequences: ", object$typeseq),
      paste0("Average length: ", lDNA),
      paste0("Number of nodes in the taxonomy: ", tot_nodes),
      paste0("Number of new nodes added to the taxonomy: ", tot_new),
      "\nTaxonomic tree and estimated parameters:",
      paste0("\t", knitr::kable(res, "simple")),
      sep = "\n")

}


