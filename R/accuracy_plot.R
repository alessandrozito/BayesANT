accuracy.plot = function(prob,correct,name="") {
  par(pty="s")
  n=length(prob)
  plot(0,type="n",xlab="cumulative prob %",ylab="cumulative correct %",xlim=c(0,100),ylim=c(0,100),las=1)
  grid()
  abline(0,1,col="gray")
  s=sort(prob,dec=F,index.return=T)
  x=cumsum(prob[s$ix])/n*100
  y=cumsum(correct[s$ix])/n*100
  lines(x,y,col="black")
  points(x[n],y[n],pch=19,cex=1,col="black")
  title(sprintf("%s\n%d items, correct %.1f %%",name,n,sum(correct)/n*100))
}

make_xy_plot <- function(prob, correct, level, id){
  s=sort(prob,dec=F,index.return=T)
  n = length(prob)
  x=cumsum(prob[s$ix])/n*100
  y=cumsum(correct[s$ix])/n*100
  df = data.frame(x, y)

  Level = paste0(id, ". ",level, " - ", as.character(round(mean(correct)*100, 1)), "%")
  l <-  ggplot2::geom_line(data = df, aes(x = x, y = y, color = Level))
  m <- ggplot2::geom_point(aes(x=x[length(x)], y=y[length(y)], color = Level))
  return(list(l, m))
}

#' @import ggplot2
plot_accuracies <- function(df_results, df_true){
  levels = colnames(df_true)[-ncol(df_true)]
  p <- ggplot2::ggplot()+
    ggplot2::xlim(0,100)+
    ggplot2::ylim(0,100)+
    ggplot2::theme_bw()+
    ggplot2::theme(aspect.ratio=1)+
    ggplot2::geom_abline(intercept = 0, slope =1, color = "grey")+
    ggplot2::xlab("cumulative prob %")+
    ggplot2::ylab("cumulative correct %") +
    ggplot2::facet_wrap(~"Model accuracies and calibration")

  for(i in 1:length(levels)){
    prob = as.numeric(df_results[, i+ length(levels)])
    correct = df_results[, i] == df_true[,i]
    f <- make_xy_plot(prob, correct, level = levels[i],id= i)
    p <- p + f[[1]] + f[[2]]
  }
  p + labs(color = "Level")
  return(p)
}

calibration_point <- function(df_results, df_true, l){
  levels = colnames(df_true)[-ncol(df_true)]
  prob = as.numeric(df_results[, l+ length(levels)])
  correct = df_results[, l] == df_true[,l]
  s=sort(prob,dec=F,index.return=T)
  n = length(prob)
  x=cumsum(prob[s$ix])/n*100
  y=cumsum(correct[s$ix])/n*100
  return(c(x[n], y[n]))
}




