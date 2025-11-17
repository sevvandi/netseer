#' Gives an error measurement for predicted graphs
#'
#' This function compares the predicted graph with the actual and comptues the
#' node and edge error as a proportion
#'
#' @param actual The ground truth or actual graph.
#' @param predicted The predicted graph.
#'
#' @return The node error and edge error as a proportion.
#'
#' @examples
#' data(syngraphs)
#' # Taking the 20th graph as the actual and the 19th graph as predicted.
#' measure_error(syngraphs[[20]], syngraphs[[19]])
#'
#'
#' @export
measure_error <- function(actual, predicted){

  node_err <- abs(igraph::vcount(actual) - igraph::vcount(predicted))/igraph::vcount(actual)
  edge_err <- abs(igraph::ecount(actual) - igraph::ecount(predicted))/igraph::ecount(actual)

                                                                                               err <- list(
    node_error = node_err,
    edge_error = edge_err
  )
  err
}
