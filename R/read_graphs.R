#' Reads graphs from a folder
#'
#' This function reads graphs from a folder to a list
#'
#' @param path_to_graphs The folder where graphs are contained.
#' @param format Formats supported by \code{igraph::read_graph}.
#'
#' @return A list of graphs in \code{igraph} format.
#'
#' @examples
#' path_to_graphs <- system.file("extdata", package = "netseer")
#' grlist <- read_graph_list(path_to_graphs = path_to_graphs, format = "gml")
#' grlist
#'
#' @export
read_graph_list <- function(path_to_graphs, format){
  if (!(format %in% c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "graphdb", "gml", "dl"))){
    stop("Invalid data format! Format needs to be one of the following: edgelist, pajek, ncol, lgl, graphml, dimacs, graphdb, gml, dl.")
  }
  # Add / to the path end if it doesn't.
  if(nchar(path_to_graphs) != "/"){
    path_to_graphs <- paste0(path_to_graphs, "/")
  }
  file_names <- list.files(path_to_graphs, pattern = paste0('*.' ,format))
  graphlist <- list()
  for(i in 1:length(file_names)){
    tempGraph <- igraph::read_graph(paste0(path_to_graphs, file_names[i]), format = format)
    graphlist[[i]] <- tempGraph
  }
  graphlist
}

