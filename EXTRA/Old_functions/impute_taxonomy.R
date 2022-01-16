imputation_via_prediction <- function(df, L, last_level, tree_branch){
  imputed_taxa = rep(NA, L)
  imputation_probs = rep(NA, L)

  for(i in 1:L){
    if(!is.na(tree_branch[i])){
      imputed_taxa[i] <- tree_branch[i]
    } else{
      if(i > last_level + 1){
        # Look at the layer below
        df = df[df[, i-1] == best, ]
      }
      # Predict
      prob <- tapply(df$predicted_probs, df[,i], sum)
      best <- names(which.max(prob))
      imputation_probs[i] <- prob[which.max(prob)]
      imputed_taxa[i] <- best
    }
  }
  cbind(data.frame(t(imputed_taxa)), data.frame(t(imputation_probs)))
}



#' Impute the taxonomy that has been removed from the original dataset
#'
#' @param x An object of class \code{BayesANT}
#' @param method Available options are \code{'prediction'} for the
#'               one-at-a-time case, and \code{'clustering'} for the many-at-a-time case (this may take several hours)
#' @return A dataset of imputed sequences
#' @export
#'
impute_taxonomy <- function(x, method = "prediction", rho=1, cores=1, adjust_weights = TRUE){

    # Number of columns in the taxonomy
  L <- ncol(x$ref_seq_missing) - 1
  if(nrow(x$ref_seq_missing)==0){
    cat("No missing labels to impute in the current taxonomy!")
    return(NULL)
  }
  if(method == "prediction"){    # <== Case 1: Prediction

    ########################################################################
    # Case 1 - impute the taxonomy via one-at-a-time case.
    ########################################################################
    doParallel::registerDoParallel(cores)
    if(x$typeseq == "aligned"){
      # Predict the aligned sequences
      out <- foreach(j = 1:nrow(x$ref_seq_missing), .combine = "rbind") %dopar% {
        impute_taxonomy_aligned(j=j,
                                ref_seq_missing = x$ref_seq_missing,
                                rho = rho,
                                M = x$M,
                                Marginalprobs = x$Marginalprobs,
                                nucl = x$nucl,
                                type_location = x$type_location,
                                adjust_weights = adjust_weights,
                                cols_to_drop = x$cols_to_drop,
                                method ="prediction")
      }
    } else if (x$typeseq == "raw"){
      # Use the k-mer decomposition
      out <- foreach(j = 1:nrow(x$ref_seq_missing), .combine = "rbind") %dopar% {
          a = impute_taxonomy_raw(j=j,
                              ref_seq_missing = x$ref_seq_missing,
                              rho = rho,
                              k= x$kmers,
                              K = x$K,
                              Marginalprobs = x$Marginalprobs,
                              nucl = x$nucl,
                              adjust_weights =  adjust_weights,
                              cols_to_drop = x$cols_to_drop,
                              method ="prediction")

      }
    }

    colnames(out) <- c(colnames(x$ref_seq_missing)[1:L], paste0("Prob_", colnames(x$ref_seq_missing)[1:L]))
    return(out)
  }
}

#a = impute_taxonomy(model, rho = 0.1, cores = 20)


