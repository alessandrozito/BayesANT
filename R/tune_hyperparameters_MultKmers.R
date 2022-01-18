

################ Method of moments tuning
Mom_dirichlet <- function(prop) {
  if (is.null(dim(prop))) {
    return(c(prop))
  } else {
    y_bar <- colMeans(prop)
    S <- mean(rowSums(prop^2))
    m <- sum(y_bar^2)
    alpha_0 <- (1 - S) / (S - m)
    if (S - m == 0) {
      return(y_bar)
    } else {
      return(alpha_0 * y_bar)
    }
  }
}


######## Get the hyperparameter counts
tune_hyperparameters_MultKmers <- function(tree, Kmer_counts, verbose = FALSE) {
  eps <- 0.01
  # Set up the hyperprior list
  hyperpriors <- vector(mode = "list", length = nrow(tree))
  cl_col <- ncol(tree)
  names(hyperpriors) <- tree[, cl_col]
  n_kmers <- ncol(Kmer_counts[[1]]$Kmers)

  ### Method 2 - Method of Moments
  props <- do.call("rbind", lapply(Kmer_counts, function(x) x[[1]]))
  props <- sweep(props + eps, 1, rowSums(props + eps), "/")

  # Note: assumption is that we have at least two levels in the hierarchy
  higher_cl <- unique(tree[, cl_col - 1])
  new_phylum <- paste0(colnames(tree)[1], "_new")

  niter <- 0
  n_clust <- length(higher_cl)
  verbose_step <- max(round(n_clust / 10), 1)

  if (verbose == T) {
    pb <- progress::progress_bar$new(
      format = "Process Taxa [:bar] :current/:total (:percent)",
      total = n_clust,
      # complete = xml2::xml_text(xml2::read_html(paste0("<x>", '&#128028;', "</x>"))),
      clear = TRUE
    )
    pb$tick(0)
  }
  tot <- 0
  for (cl in higher_cl) {
    niter <- niter + 1
    if (grepl("_new$", cl)) {
      if (cl == new_phylum) {
        clusters <- tree[, cl_col]
      } else {
        # Find the parent level
        parent_level <- which(tree[which(tree[, cl_col - 1] == cl), ] == cl)[1] - 1
        parent <- tree[, parent_level][tree[, cl_col - 1] == cl]
        clusters <- tree[, cl_col][tree[, parent_level] == parent]
      }
      select <- (names(Kmer_counts) %in% clusters)
      indexes <- which(names(hyperpriors) == cl)
    } else {
      # Cluster not new
      clusters <- tree[, cl_col][tree[, cl_col - 1] == cl]
      indexes <- which(names(hyperpriors) %in% clusters)
      select <- (names(Kmer_counts) %in% clusters)
    }

    alpha_mat <- t(as.matrix(Mom_dirichlet(props[select, ])))

    for (i in indexes) {
      hyperpriors[[i]] <- unname(alpha_mat)
    }

    if (niter %% verbose_step == 0 & verbose == TRUE) {
      pb$tick(verbose_step)
      tot <- tot + verbose_step
    }
  }
  if (verbose & length(higher_cl) - tot > 0) {
    pb$tick(length(higher_cl) - tot)
  }
  return(hyperpriors)
}
