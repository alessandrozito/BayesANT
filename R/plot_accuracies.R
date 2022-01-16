#' Calibration plot for the output of a BayesANT prediction
#'
#' @param prediction An object of class \code{data.frame} that is the output of the function \code{predict.BayesANT}
#' @param test_data An object of class \code{c('data.frame','BayesANT.data')} that contains the true taxonomic annotations.
#'                  It is useful to process this first via the function \code{add_novelty_to_test_data}
#'
#' @import ggplot2
#' @export
plot_accuracies <- function(prediction, test_data) {
  levels <- colnames(test_data)[-ncol(test_data)]
  p <- ggplot2::ggplot() +
    ggplot2::xlim(0, 100) +
    ggplot2::ylim(0, 100) +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey") +
    ggplot2::xlab("cumulative prob %") +
    ggplot2::ylab("cumulative correct %") +
    ggplot2::facet_wrap(~"Model accuracies and calibration")

  for (i in 1:length(levels)) {
    prob <- as.numeric(prediction[, i + length(levels)])
    correct <- prediction[, i] == test_data[, i]
    f <- make_xy_plot(prob, correct, level = levels[i], id = i)
    p <- p + f[[1]] + f[[2]]
  }
  p + labs(color = "Level")
  return(p)
}


make_xy_plot <- function(prob, correct, level, id) {
  s <- sort(prob, decreasing = F, index.return = T)
  n <- length(prob)
  x <- cumsum(prob[s$ix]) / n * 100
  y <- cumsum(correct[s$ix]) / n * 100
  df <- data.frame(x, y)

  Level <- paste0(id, ". ", level, " - ", as.character(round(mean(correct) * 100, 1)), "%")
  l <- ggplot2::geom_line(data = df, aes(x = x, y = y, color = Level))
  m <- ggplot2::geom_point(aes(x = x[length(x)], y = y[length(y)], color = Level))
  return(list(l, m))
}
