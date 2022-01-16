samplePY <- function(size, alpha, sigma) {
  if(size == 1){
    return(1)
  }
  # Sample a Pitman-Yor of a given size
  counts <- c(1)
  for (n in 2:size) {
    K <- length(counts)
    prob_new <- (alpha + K * sigma) / (n - 1 + alpha)
    prob_old <- (counts - sigma) / (n - 1 + alpha)

    # Sample the observation
    x <- sample(x = c(0, 1:K), size = 1, prob = c(prob_new, prob_old), replace = FALSE)

    if (x == 0) {
      counts <- c(counts, 1)
    } else {
      counts[x] <- counts[x] + 1
    }
  }
  return(counts)
}


#sampleDM <- function(size, sigma, H) {
#  if(sigma>0){
#    cat("sigma must be negative!")
#  } else {
#    H = as.integer(H)
#    alpha = -H*sigma
#    samplePY(size, alpha = alpha,sigma = sigma)
#  }
#}

