plot_similarity <- function(DNA, seqtype = "aligned", kmers = 5, min_score = 0, taxa_to_plot = NULL, colors = c("#2960F5", "#BDEAF5")){
  dnabin = ape::as.DNAbin(x = stringr::str_split(DNA, ""))
  if(seqtype == "aligned"){
    M =  1 - ape::dist.dna(x = dnabin, model = "raw", as.matrix = TRUE)
  } else if (seqtype == "not aligned"){
    M = 1 - as.matrix(kdistance(x = dnabin, k = kmers))
    rownames(M) <- colnames(M) <- 1:ncol(M)
    #M<- t(apply(M, 2, rev))

  }
  longData<-reshape2::melt(M[,rev(1:ncol(M))])

  p <- ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    #scale_fill_gradient2(low="yellow",midpoint = 0.5, mid ="orange",
    #                    high="red3", limits = c(0,1)) +
    scale_y_reverse()+
    scale_fill_distiller("Similarity", palette = "YlOrRd",limits = c(min_score,1), direction = 1)+
    theme_minimal() +
    theme(aspect.ratio = 1,
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    #scale_y_reverse() +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 15, title = "")) +
    theme( legend.margin=margin(0,0,0,0),
           legend.box.margin=margin(-10,-10,-10,-10))

  if(!is.null(taxa_to_plot)){
    ## Add the rectangle with the taxonomy. One for each column
    for(j in 1:ncol(taxa_to_plot)){
      borders <- which(!duplicated(taxa_to_plot[,j])) - 0.5
      p <- p + annotate(geom='rect', xmin= borders[-length(borders)],
                        ymin= borders[-length(borders)],
               xmax=borders[-1], ymax=  borders[-1],
               fill="transparent", col=colors[j], lwd=0.65)
    }
  }
  return(p)
}

#YlOrRd
