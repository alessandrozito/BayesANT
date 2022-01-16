#' Relabel the true test annotations based on a reference library used to train a BayesANT model
#'
#' @param data An object of class \code{c('data.frame', 'BayesANT.data')}
#' @param test_data An object of class \code{c('data.frame', 'BayesANT.data')}
#'
#' @return An object of class \code{c('data.frame', 'BayesANT.data')}
#' @export
add_novelty_to_test_data <- function(data, test_data) {
  if (!all(class(data) == c("data.frame", "BayesANT.data"))) {
    stop("data must be of class 'data.frame', 'BayesANT.data'")
  }
  if (!all(class(test_data) == c("data.frame", "BayesANT.data"))) {
    stop("test_data must be of class 'data.frame', 'BayesANT.data'")
  }

  df_true <- test_data

  # Get the names of the levels
  level_names <- names(data)[-ncol(data)]
  max_levels <- length(level_names)
  for (l in 1:max_levels) {
    # Find the indexes of new species
    id <- !(df_true[, l] %in% unique(data[, l]))
    n_replace <- sum(id)
    if (n_replace > 0) {
      if (l == 1) {
        name_clust <- paste0(level_names[l], "_new")
        df_true[id, l:max_levels] <- matrix(rep(name_clust, n_replace * max_levels),
          nrow = n_replace, ncol = max_levels
        )
      } else {
        for (i in which(id == TRUE)) {
          if (!grepl("_new", df_true[i, l - 1])) {
            name_clust <- paste0(df_true[i, l - 1], "_", level_names[l], "_new")
            df_true[i, l:max_levels] <- rep(name_clust, max_levels - l + 1)
          }
        }
      }
    }
  }
  class(df_true) <- c("data.frame", "BayesANT.data")
  return(df_true)
}
