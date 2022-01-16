##### Method of moments tuning

#' Extract the proportions of the given Nucleotide counts
#'
#' @param Nucl_counts nulcelotide counts list
#' @param eps small tweak for the zero counts
#'
#' @return
get_proportions <- function(Nucl_counts, eps = 0, verbose = FALSE){

  Proportions <- vector("list", length = ncol(Nucl_counts[[1]]))

  verbose_step <- round(length(Proportions)/10, 0)
  if(verbose){
  pb <- progress::progress_bar$new(format = "Process locations [:bar] :current/:total (:percent)",
                                   total = length(Proportions),
                                   #complete = xml2::xml_text(xml2::read_html(paste0("<x>", '&#128028;', "</x>"))),
                                   clear = TRUE)

    pb$tick(0)
  }

  tot <- 0
  for(i in 1:length(Proportions)){

    #if(i%%verbose_step==0 & verbose == TRUE){
    #  cat(paste0("Number of locations processed: ", i, "\n"))
    #}

    Proportions[[i]] <- t(sapply(Nucl_counts, function(x) (x[,i] + eps)/(sum(x[,i] + eps)), simplify = "matrix"))
    if(i%%verbose_step==0 & verbose == TRUE){
      pb$tick(verbose_step)
      tot <- tot + verbose_step
    }
  }
  if(verbose){
    pb$tick(length(Proportions) - tot)
  }
  return(Proportions)
}

################ Method of moments tuning
Mom_dirichlet <- function(prop){
  if(is.null(dim(prop))){
    return(c(prop))
  } else {
    y_bar <- colMeans(prop)
    S <- mean(rowSums(prop^2))
    m <- sum(y_bar^2)
    alpha_0 <- (1-S)/(S-m)
    if(S-m == 0){
      return(y_bar)
    } else {
      return(alpha_0*y_bar)
    }
  }
}


######## Get the hyperparameter counts
tune_hyperparameters <- function(tree, Nucl_counts, mixture_weight = NULL, verbose = FALSE){

  # Set up the hyperprior list
  hyperpriors = vector(mode = "list", length = nrow(tree))
  cl_col = ncol(tree)
  names(hyperpriors) = tree[, cl_col]
  length_seq <- ncol(Nucl_counts[[1]])
  basis = nrow(Nucl_counts[[1]])

  if(is.null(mixture_weight)){mixture_weight = 1/basis}

  ### Method 2 - Method of Moments

    ## Get the proportions
    eps <- 0.01
    Proportions <- get_proportions(Nucl_counts, eps = eps, verbose = verbose)

    # Note: assumption is that we have at least two levels in the hierarchy
    higher_cl = unique(tree[, cl_col - 1])
    new_phylum =paste0(colnames(tree)[1], "_new")

    niter = 0
    n_clust = length(higher_cl)
    verbose_step = round(n_clust/10)
    if(verbose == T){

    pb <- progress::progress_bar$new(format = "Process Taxa [:bar] :current/:total (:percent)",
                                     total = n_clust,
                                     #complete = xml2::xml_text(xml2::read_html(paste0("<x>", '&#128028;', "</x>"))),
                                     clear = TRUE)
    pb$tick(0)
    }
    tot <- 0
    for(cl in higher_cl){
      niter = niter + 1
      #if(niter %% verbose_step == 0 & verbose == TRUE){
      #  cat(paste0("Number of taxa processed: ", niter, " out of ", n_clust, "\n"))
      #}
      if(grepl("_new$", cl)){

        if(cl == new_phylum){
          clusters = tree[, cl_col]
        } else {
          # Find the parent level
          parent_level <- which(tree[which(tree[,cl_col-1] == cl), ] == cl)[1] - 1
          parent <- tree[, parent_level][tree[, cl_col - 1] == cl]
          clusters = tree[, cl_col][tree[, parent_level] == parent]
        }
        select = (rownames(Proportions[[1]]) %in% clusters)
        indexes = which(names(hyperpriors) == cl)

      } else {
        # Families
        clusters = tree[, cl_col][tree[, cl_col-1] == cl]
        indexes = which(names(hyperpriors) %in% clusters)
        select = (rownames(Proportions[[1]]) %in% clusters)
      }

      # Aggregate
      alpha_mat = t(do.call("rbind", lapply(Proportions, function(x) Mom_dirichlet(x[select, ]))))
      for(i in indexes){hyperpriors[[i]] = unname(alpha_mat)}

      if(niter %% verbose_step == 0 & verbose == TRUE){
        pb$tick(verbose_step)
        tot <- tot + verbose_step
      }

    }
    if(verbose==T){
      pb$tick(length(higher_cl) - tot)
    }
  return(hyperpriors)
}





