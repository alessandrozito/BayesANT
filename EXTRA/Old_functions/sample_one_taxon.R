sample_one_taxon <- function(Priorprobs){
  depth = (ncol(Priorprobs) - 2)/2
  prediction_probs <- rep(NA, depth)
  prediction_names <- rep(NA, depth)

  for(i in 1:length(prediction_probs)){
    if(i == 1){
      df = Priorprobs
    } else {
      df = Priorprobs[Priorprobs[, i-1] == sampled_taxon, ]
    }
    prob <- tapply(df$predicted_probs, df[,i], sum)
    sampled_taxon <- sample(x = names(prob), size = 1, prob = prob)  # Sample instead of choosing the max
    prediction_probs[i] <- prob[which(names(prob)==sampled_taxon)]
    prediction_names[i] <- sampled_taxon
  }
  out <- c(unname(prediction_names), unname(prediction_probs))
  return(out)
}
