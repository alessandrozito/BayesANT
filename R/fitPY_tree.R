##### Empirical Bayes estimates for the Pitman-Yor process

logEPPF_PY = function(alpha, sigma, frequencies) {
  if (any(sigma < 0, alpha <= -sigma + 1e-04)) {
    return(-Inf)
  }

  # Sample size
  n = sum(frequencies)
  # Number of frequencies
  K = length(frequencies)
  # Loglikelihood
  loglik = sum(log(alpha + 1:(K - 1) * sigma)) - lgamma(alpha + n) + lgamma(alpha + 1) + sum(lgamma(frequencies - sigma)) - K * lgamma(1 - sigma)

  loglik
}

logEPPF_PY_pooled = function(alpha, sigma, freq_list){
  # Pooled Loglikelihood (all freq vectors in freq_list are IID from the same PY)
  loglik_pooled = sum(unlist(lapply(freq_list, FUN = function(f) logEPPF_PY(alpha, sigma, f))))
  loglik_pooled
}

max_EPPF_PY = function(frequencies) {
  start = c(1, 0.5) # Initialization of the maximization algorithm
  out = nlminb(
    start = start,
    function(param) -logEPPF_PY(alpha = param[1], sigma = param[2], frequencies = frequencies),
    lower = c(-Inf, 1e-16), upper = c(Inf, 1 - 1e-10)
  )
  return(out)
}

max_EPPF_PY_pooled = function(freq_list){
  start = c(1, 0.5) # Initialization of the maximization algorithm
  out = nlminb(
    start = start,
    function(param) -logEPPF_PY_pooled(alpha = param[1], sigma = param[2], freq_list = freq_list),
    lower = c(-Inf, 1e-16), upper = c(Inf, 1 - 1e-10)
  )
  return(out)
}

fitPY <- function(freq) {
  if (is.list(freq)) {
    # If freq is a list, fit the Pooled model
    fit = max_EPPF_PY_pooled(freq)
    out = list(frequencies = freq, param = fit$par, loglik = -fit$objective)
    return(out)
  } else {
    # If freq is a vector, fit a PY for a single vector of counts
    fit = max_EPPF_PY(frequencies = freq)
    out = list(frequencies = freq, param = fit$par, loglik = -fit$objective)
    return(out)
  }
}

#' Fit the Pitman-Yor parameters
#'
#' @param tree A tree that comes out of the simulate_tree
#'
#' @return
fitPY_tree <- function(tree){
  # Assume that the last column contains the sequences
  PYpars = vector(mode = "list", length = ncol(tree)-1)
  names(PYpars) = colnames(tree)[-ncol(tree)]

  # Level 1
  PYpars[[1]] =  fitPY(table(tree[,1]))

  # Other levels
  for(l in 2:length(PYpars)){
    tab = table(tree[, l-1], tree[, l])
    freq_list = lapply(1:nrow(tab), function(i) tab[i, ][tab[i, ]>0])
    names(freq_list) = rownames(tab)
    PYpars[[l]] = fitPY(freq_list)
  }
  return(PYpars)
}


