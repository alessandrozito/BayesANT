obtain_prediction <- function(Priorprobs) {
  depth <- (ncol(Priorprobs) - 2) / 2
  prediction <- rep(NA, depth)
  prediction_names <- rep(NA, depth)

  for (i in 1:length(prediction)) {
    if (i == 1) {
      df <- Priorprobs
    } else {
      df <- Priorprobs[Priorprobs[, i - 1] == best, ]
    }
    prob <- tapply(df$predicted_probs, df[, i], sum)
    best <- names(which.max(prob))
    prediction[i] <- prob[which.max(prob)]
    prediction_names[i] <- best
  }
  out <- c(unname(prediction_names), unname(prediction))
  return(out)
}
