get_nucleotide_counts <- function(nodes, DNAseq, nucl, verbose = F) {

  # Extract the unique families
  unique_nodes <- sort(unique(nodes))

  # Initialize an empty list
  Nucl_counts <- vector(mode = "list", length = length(unique_nodes))
  names(Nucl_counts) <- unique_nodes

  # Progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Nucleotide counts [:bar] :current/:total (:percent)",
      total = length(unique_nodes),
      # complete = xml2::xml_text(xml2::read_html(paste0("<x>", '&#128028;', "</x>"))),
      clear = TRUE
    )

    pb$tick(0)
  }
  pb_step <- round(length(unique_nodes) / 100)
  tot <- 0
  for (i in 1:length(unique_nodes)) {
    indexes <- which(unique_nodes[i] == nodes)

    if (length(indexes) > 1) {
      matDNA <- do.call(rbind, stringr::str_split(DNAseq[indexes], ""))
      matDNA[!matDNA %in% nucl] <- "-"
    } else {
      matDNA <- t(as.matrix(stringr::str_split(DNAseq[indexes], "")[[1]]))
      matDNA[!matDNA %in% nucl] <- "-"
    }

    Nucl_counts[[i]] <- build_CountMatrix(matDNA, nucl = nucl)
    # Print the progress bar
    if (i %% pb_step == 0 & verbose) {
      pb$tick(pb_step)
      tot <- tot + pb_step
    }
  }
  if (verbose) {
    pb$tick(length(unique_nodes) - tot)
  }

  return(Nucl_counts)
}
