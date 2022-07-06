#' Main function to call the BayesANT (BAYESiAn Nonparametric Taxonomic)
#' classifier. The function construct an object of class \code{'BayesANT'},
#' which can be used for predictions.
#'
#' @param data An object of class \code{c('data.frame', 'BayesANT.data')}
#'             containing the training library. Needs to be loaded with the
#'             function \code{read.BayesANT.data}.
#' @param typeseq Format of sequences used to train the classifier.
#'                Options are \code{'not aligned'} and \code{'aligned'}.
#'                Default is \code{typeseq = 'aligned'}. In this last case,
#'                the function returns an error if the sequences are not of the
#'                same length.
#' @param kmers Length of a substring under a k-mer decomposition.
#'              Valid only if \code{typeseq = 'not aligned'}.
#'              Default is  \code{kmers = 5}, which is also the recommended
#'              choice. The maximum choice allowed is \code{kmers = 8}
#' @param type_location How to model the loci in an aligned set of DNA
#'                      sequences. Valid only if \code{typeseq = 'aligned'}.
#'                      Options are \code{'single'}, which refers to the
#'                      single-location multinomial kernel, and  \code{'pairs'},
#'                      which is the 2-mer multinomial kernel.  Default is
#'                      \code{type_location = 'single'}
#' @param usegap Whether to include the alignment gap \code{"-"} among the
#'               nucleotides in the aligned multinomial kernel. Default is
#'               \code{usegap = FALSE}
#' @param newtaxa Whether to account for new taxa when constructing the
#'                classifier. Default is \code{newtaxa = TRUE}. If
#'                \code{newtaxa = FALSE},no potential unobserved branches are
#'                included in the taxonomy.
#' @param save_nucl Save the nucleotide counts and the hyperparameters of
#'                  the model in a list. Default is \code{save_nucl = FALSE}.
#'                  Setting it to \code{TRUE} might be heavy to store, so use
#'                  only if strictly needed.
#' @param adjust_Kmer_length Adjust the likelihood contribution of a k-mer
#'                           decomposition based on its length.
#'                           Valid for the option \code{typeseq ='not aligned'}.
#'                           Default is \code{adjust_Kmer_length = TRUE}.
#' @param verbose Monitor the steps adopted to train the algorithm.
#'                Default is \code{verbose = TRUE}.
#'
#' @return An object of class \code{BayesANT}. We return a \code{list}
#'         containing the following quantities:
#'         \describe{
#'          \item{data}{Dataset used for training.},
#'          \item{data_missing}{Dataset containing sequences with missing
#'             annotations.}
#'          \item{missing_taxa}{Indeces of the sequences with missing values.}
#'          \item{typeseq}{Type of sequences used to train the classifier.}
#'          \item{type_location}{Type of location to train the classifier.}
#'          \item{newtaxa}{Whether new taxa are included in the classification.}
#'          \item{level_names}{Names of taxonomic ranks}
#'          \item{nucl}{Nucleotides detected.}
#'          \item{kmers}{Number of kmers selected to build the classifier.}
#'          \item{ParameterMatrix}{Matrix that stores model parameters.}
#'          \item{Nucl_counts}{List of counts of the nucleotides at every leaf.}
#'          \item{PYpars}{Estimated Pitman-yor parameters.}
#'          \item{Priorprobs}{Prior probabilities selected for the model.}
#'          \item{hyperparameters}{List of model hyperparameters.}
#'          \item{leaves}{Name of the taxonomic leaves.}
#'          \item{adjust_Kmer_length}{Wheter to adjust for sequence length.}
#'          \item{sequences_length}{Length of each sequence in the data.}
#'         }
#' @export
BayesANT <- function(data,
                     typeseq = "aligned",
                     type_location = "single",
                     kmers = 5,
                     newtaxa = TRUE,
                     usegap = FALSE,
                     adjust_Kmer_length = TRUE,
                     save_nucl = FALSE,
                     verbose = TRUE) {
  # Check if parameters of the function as correctly specified.
  stopifnot(
    class(data) == c("data.frame", "BayesANT.data"),
    typeseq %in% c("aligned", "not aligned"),
    type_location %in% c("single", "pairs"),
    kmers > 0 & kmers <= 8,
    is.logical(newtaxa),
    is.logical(usegap),
    is.logical(save_nucl),
    is.logical(verbose)
  )

  if (verbose) {
    cat(
      "Welcome to BayesANT",
      rep(xml2::xml_text(
        xml2::read_html(paste0("<x>", "&#128028;", "</x>"))
      ), 4), " \n"
    )
  }



  # Save the Identifiers of the sequences.
  IDs <- rownames(data)

  # Exclude the observations with missing values.
  # We store the index of the incomplete cases if one is interested
  # in clustering
  missing_taxa <- !stats::complete.cases(data)
  data_missing <- data[missing_taxa, ]
  # Restrict to the fully annotated sequences
  data <- data[!missing_taxa, ]
  if (verbose) {
    if (sum(missing_taxa) == 0) {
      cat("No missing annotations detected. \n")
    } else {
      cat(paste0(
        "Removing ", sum(missing_taxa),
        " observations from the taxonomy due to missingness. \n"
      ))
    }
  }

  # Extract the information
  lowest_level <- ncol(data) - 1 ## The last column in the data DNA sequences
  nodes <- data[, lowest_level]
  DNAseq <- data[, lowest_level + 1]

  if (usegap) {
    nucl <- c("-", "A", "C", "G", "T")
  } else {
    # Do not use the '-' among the nucleotides
    nucl <- c("A", "C", "G", "T")
  }

  # Do one check on the length of the DNA sequences if the
  # specified kernel is aligned
  seqs_length <- unique(stringr::str_length(DNAseq))
  if ((length(seqs_length) != 1) & typeseq == "aligned") {
    stop("DNA sequences of different lengths are not allowed when
         typeseq = 'aligned'")
  }

  #--------------------------------------------------------------
  # calculate the Pitman-Yor weights at every node.
  #--------------------------------------------------------------
  if (verbose) {
    cat("Estimating the Pitman-Yor parameters: Now. \n")
  }
  PYpars <- suppressWarnings(fitPY_tree(data))
  if (verbose) {
    cat("Estimating the Pitman-Yor parameters: Done. \n")
  }

  #--------------------------------------------------------------
  # Extract the prior probabilities given by the Pitman-Yor processes
  #--------------------------------------------------------------
  if (verbose) {
    cat("Computing prior probabilities: Now \n")
  }

  if (newtaxa) {
    Priorprobs <- get_Priorprobs_PY(data, PYpars)
  } else {
    Priorprobs <- get_Priorprobs_noNewSpecies(data, PYpars)
  }

  leaves <- Priorprobs[, lowest_level]
  if (verbose) {
    cat("Computing prior probabilities: Done. \n")
  }

  #--------------------------------------------------------------
  # Sequence processing and parameter estimation
  #--------------------------------------------------------------

  #------------------------------------------------- Aligned sequences
  if (typeseq == "aligned") {
    if (type_location == "single") {
      # Use the classic multinomial kernel
      if (verbose) {
        cat("Extracting nucleotide location counts from sequences: Now. \n")
      }
      Nucl_counts <- get_nucleotide_counts(nodes, DNAseq, nucl)
      if (verbose) {
        cat("Extracting nucleotide location counts from sequences: Done. \n")
      }
    } else if (type_location == "pairs") {
      # Pairs of nucleotides under the classic multinomial kernel
      nucl <- paste0(rep(nucl, each = length(nucl)), rep(nucl, length(nucl)))
      if (verbose) {
        cat("Extracting 2mer location counts from sequences: Now. \n")
      }
      Nucl_counts <- get_2mer_counts(nodes, DNAseq, nucl)
      if (verbose) {
        cat("Extracting 2mer location counts from sequences: Done. \n")
      }
    }

    #------------------------------------------------- Not aligned sequences
  } else if (typeseq == "not aligned") {
    # The sequences are not aligned. We need to use the kmer decomposition
    if (verbose) {
      cat("Extracting kmers from sequences: Now. \n")
    }
    Nucl_counts <- get_Kmer_counts(
      k = kmers, nodes, DNAseq, nucl,
      adjust_Kmer_length = adjust_Kmer_length
    )
    if (verbose) {
      cat("Extracting kmers from sequences: Done. \n")
    }
  }

  #--------------------------------------------------------------
  # Calculating hyperparameters via method of the moments.
  #--------------------------------------------------------------
  if (verbose) {
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

  if (verbose) {
    cat("Tuning hyperparameters: Done. \n")
  }

  #--------------------------------------------------------------
  # Aggregate all parameters
  #--------------------------------------------------------------
  if (verbose) {
    cat("Aggregating parameters: Now \n")
  }
  if (typeseq == "aligned") {
    ParameterMatrix <- build_ParameterMatrix(
      leaves,
      Nucl_counts,
      hyperparameters
    )
  } else if (typeseq == "not aligned") {
    ## Build the matrix from the kmers
    ParameterMatrix <- build_ParameterMatrix_from_Kmers(
      leaves,
      Nucl_counts,
      hyperparameters
    )
  }
  if (verbose) {
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
