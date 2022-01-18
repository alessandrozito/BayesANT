#' Call to the BayesANT Classifier.
#'
#' @param data An object of class \code{c('data.frame', 'BayesANT.data')}. Needs to be loaded with the function \code{read.BayesANT.data}.
#' @param typeseq Type of sequences used to train the classifier. Options are \code{'not aligned'} and \code{'aligned'}. Default is \code{typeseq = 'aligned'}
#' @param kmers How many kmers to take for the k-mer decomposition. Valid only if \code{typeseq = 'not aligned'}. Default is  \code{kmers = 5}
#' @param type_location How to model the loci in an aligned set of DNA sequences. Valid only if \code{typeseq = 'aligned'}.
#'                      Options are \code{'single'}, which refers to the single-location multinomial kernel, and  \code{'pairs'}, which is the 2-Product multinomial kernel. Default is \code{type_location = 'single'}
#' @param usegap Whether to include the alignment gap \code{"-"} among the nucleotides in the aligned multinomial kernel. Default is \code{usegap = FALSE}
#' @param newtaxa Include new taxa when constructing the classifier. Default is \code{newtaxa = TRUE}
#' @param save_nucl Save the nucleotide counts and the hyperparameters of the model in a list. Default is \code{save_nucl = FALSE}. Setting it to \code{TRUE} might severely increase the storage requirements, so be careful.
#' @param adjust_Kmer_length Adjust the likelihood contribution of a kmer decomposition based on its length. Valid for the option  \code{typeseq ='not aligned'}.
#' @param verbose Monitor the state of the algorithm. Default is \code{verbose = TRUE}
#'
#' @return An object of class \code{BayesANT}
#' @export
#'
BayesANT <- function(data,
                     typeseq = "aligned",
                     type_location = "single",
                     kmers = 5,
                     newtaxa = TRUE,
                     usegap = FALSE,
                     adjust_Kmer_length = TRUE,
                     save_nucl = FALSE,
                     verbose = TRUE) {
  if (verbose == TRUE) {
    cat("Welcome to BayesANT", rep(xml2::xml_text(xml2::read_html(paste0("<x>", "&#128028;", "</x>"))), 4), " \n")
  }

  ## Check if the reference library is an object of class "BayesANT.data"
  if (!all(class(data) == c("data.frame", "BayesANT.data"))) {
    stop("data must be of class 'data.frame', 'BayesANT.data'")
  }

  IDs <- rownames(data)

  # Exclude the observations with missing values.
  # We store the index of the incomplete cases if one is interested in clustering
  missing_taxa <- !stats::complete.cases(data)
  data_missing <- data[missing_taxa, ]
  # Restrict to the fully annotated sequences
  data <- data[!missing_taxa, ]
  if (verbose == TRUE) {
    if (sum(missing_taxa) == 0) {
      cat("No missing annotations detected. \n")
    } else {
      cat(paste0("Removing ", sum(missing_taxa), " observations from the taxonomy due to missingness. \n"))
    }
  }

  # Extract the information
  lowest_level <- ncol(data) - 1 ## The last column in the data are the DNA sequences
  nodes <- data[, lowest_level]
  DNAseq <- data[, lowest_level + 1]

  if (usegap) {
    nucl <- c("-", "A", "C", "G", "T")
  } else {
    # Do not use the '-' among the nucleotides
    nucl <- c("A", "C", "G", "T")
  }

  # Do one check on the length of the DNA sequences if the specified kernel is aligned
  seqs_length <- unique(stringr::str_length(DNAseq))
  if ((length(seqs_length) != 1) & typeseq == "aligned") {
    stop("DNA sequences of different lengths are not allowed when typeseq = 'aligned'")
  }

  ##### STEP 1 - extract the Pitman-Yor weights
  if (verbose == TRUE) {
    cat("Estimating the Pitman-Yor parameters: Now. \n")
  }
  PYpars <- suppressWarnings(fitPY_tree(data))
  if (verbose == TRUE) {
    cat("Estimating the Pitman-Yor parameters: Done. \n")
  }

  ##### STEP 2 - extract the Prior probabilities given by the Pitman-Yor processes
  if (verbose == TRUE) {
    cat("Computing prior probabilities: Now \n")
  }

  if (newtaxa == TRUE) {
    Priorprobs <- get_Priorprobs_PY(data, PYpars)
  } else {
    Priorprobs <- get_Priorprobs_noNewSpecies(data, PYpars)
  }

  leaves <- Priorprobs[, lowest_level]
  if (verbose == TRUE) {
    cat("Computing prior probabilities: Done. \n")
  }

  ##### STEP 3 - Processing the sequences
  # Process the nucleotides. In particular, we can either use the
  # kmer decomposition or look at the individual locations

  if (typeseq == "aligned") {
    if (type_location == "single") {
      ## Use the classic multinomial kernel
      if (verbose == TRUE) {
        cat("Extracting nucleotide location counts from sequences: Now. \n")
      }
      Nucl_counts <- get_nucleotide_counts(nodes, DNAseq, nucl)
      if (verbose == TRUE) {
        cat("Extracting nucleotide location counts from sequences: Done. \n")
      }
    } else if (type_location == "pairs") {
      # Pairs of nucloetides under the classic multinomial kernel
      prm <- gtools::permutations(n = length(nucl), r = 2, v = nucl, repeats.allowed = TRUE)
      nucl <- apply(prm, 1, function(x) paste0(x, collapse = ""))
      if (verbose == TRUE) {
        cat("Extracting 2mer location counts from sequences: Now. \n")
      }
      Nucl_counts <- get_2mer_counts(nodes, DNAseq, nucl)
      if (verbose == TRUE) {
        cat("Extracting 2mer location counts from sequences: Done. \n")
      }
    } else {
      stop(paste0("typeseq 'aligned' can only work with location = 'single' or 'pairs'"))
    }
  } else if (typeseq == "not aligned") {
    # The sequences are not aligned. We need to use the kmer decomposition
    if (verbose == TRUE) {
      cat("Extracting kmers from sequences: Now. \n")
    }
    Nucl_counts <- get_Kmer_counts(
      k = kmers, nodes, DNAseq, nucl,
      adjust_Kmer_length = adjust_Kmer_length
    )
    if (verbose == TRUE) {
      cat("Extracting kmers from sequences: Done. \n")
    }
  } else {
    stop(cat("Invalid typeseq specified. Options allowed are 'aligned' and 'not aligned' \n"))
  }

  ##### STEP 4 - computing the hyperparameters
  # Use the method of the moments to extract the hyperpparameters
  if (verbose == TRUE) {
    cat("Tuning hyperparameters: Now. \n")
  }

  if (typeseq == "aligned") {
    hyperparameters <- tune_hyperparameters(
      tree = Priorprobs[, 1:lowest_level],
      Nucl_counts = Nucl_counts, verbose = verbose
    )

  } else if (typeseq == "not aligned") {
    hyperparameters <- tune_hyperparameters_MultKmers(
      tree = Priorprobs[, 1:lowest_level],
      Kmer_counts = Nucl_counts,
      verbose = verbose
    )
  }

  if (verbose == TRUE) {
    cat("Tuning hyperparameters: Done. \n")
  }


  ##### STEP 4 - Aggregating the parameters
  if (verbose == TRUE) {
    cat("Aggregating parameters: Now \n")
  }
  if (typeseq == "aligned") {
    ParameterMatrix <- build_ParameterMatrix(leaves, Nucl_counts, hyperparameters)
  } else if (typeseq == "not aligned") {
    ## Build the matrix from the kmers
    ParameterMatrix <- build_ParameterMatrix_from_Kmers(leaves, Nucl_counts, hyperparameters)
  }
  if (verbose == TRUE) {
    cat("Aggregating parameters: Done \n")
  }

  if (save_nucl == FALSE) {
    # eliminate the hyperpatameters and the nucleotide counts to save space
    Nucl_counts <- NULL
    hyperparameters <- NULL
  }

  out <- list(
    "data" = data,
    "data_missing" = data_missing,
    "missing_taxa" = missing_taxa,
    "typeseq" = typeseq,
    "type_location" = type_location,
    "newtaxa" = newtaxa,
    "level_names" = colnames(data)[-ncol(data)],
    "nucl" = nucl,
    "kmers" = kmers,
    "ParameterMatrix" = ParameterMatrix,
    "Nucl_counts" = Nucl_counts,
    "PYpars" = PYpars,
    "Priorprobs" = Priorprobs,
    "hyperparameters" = hyperparameters,
    "leaves" = leaves,
    "adjust_Kmer_length" = adjust_Kmer_length,
    "sequences_length" = seqs_length
  )

  class(out) <- "BayesANT"
  return(out)
}
