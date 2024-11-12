#' Predicts a graph from a time series of graphs.
#'
#' This function predicts the graph at a future time step using a time series of
#' graphs.
#'
#' @param graphlist A list of graphs in igraph format.
#' @param formulation Formulation 2 includes an additional condition constraining total
#' edges by the predicted value. Formulation 1 does not have that constraint. Formulation 2
#' gives more realistic graphs due to that constraint. Default is set to \code{2}.
#' @param conf_level1 A value between 50 and 100 denoting the confidence interval
#' for the number of predicted nodes in the graph. If set to \code{NULL} the predicted
#' graph has the mean number of predicted nodes. If set to \code{80} for example,
#' there would be 3 predicted graphs. One with mean number of predicted nodes, and
#' the other two with the number of nodes corresponding to lower and upper
#' confidence bounds.
#' @param conf_level2 The upper confidence bound for the degree distribution. Default
#' set to \code{90}.
#' @param dense_opt If set to \code{2} the dense option in R package \code{lpSolve}
#' will be used.
#' @param weights_opt Weights option ranging from 1 to 6 used for different edge weight
#' schemes. Weights option 1 uses uniform weights for all edges. Option 2 uses binary
#' weights. If the edge existed in a past graph, then weight is set to 1. Else set to
#' 0. All possible new edges are assigned weight 1. Option 3 is a more selective
#' version. Option 4 uses proportional weights according to the history. Option 5 uses
#' proportional weights, but as the network is more in the past, it gives less weight.
#' Option 5 uses linearly decaying proportional weights. Option 6 uses harmonically decaying
#' weights. That is the network at \code{T} is given weight 1,  \code{T-1}
#' is given weight 1/2 and so on. Default is set to \code{6}.
#' @param weights_param The weight given for possible edges from new vertices. Default
#' set to \code{0.001}.
#' @param h The prediction time step. Default is \code{ h = 1}.
#'
#'
#' @return A list of predicted graphs. If \code{conf_level1} is not \code{NULL}, then
#' 3 graphs are returned one with the mean number of predicted nodes and the other 2
#' with the number of nodes equal to the lower and upper bound values of prediction.
#' If If \code{conf_level1} is \code{NULL}, only the mean predicted graph is returned.
#'
#' @examples
#' set.seed(2024)
#' edge_increase_val <- new_nodes_val <- del_edge_val <- 0.1
#' graphlist <- list()
#' graphlist[[1]] <- gr <-  igraph::sample_pa(5, directed = FALSE)
#' for(i in 2:15){
#'   gr <-  generate_graph_exp(gr,
#'                         del_edge = del_edge_val,
#'                         new_nodes = new_nodes_val,
#'                         edge_increase = edge_increase_val )
#'   graphlist[[i]] <- gr
#' }
#' grpred <- predict_graph(graphlist[1:15], conf_level2 = 90, weights_opt = 6)
#' grpred
#'  \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   if (!inherits(future::plan(), "sequential")) future::plan(sequential)
#'   }
#'
#' @importFrom dplyr pull summarize mutate group_by n full_join filter arrange
#' @importFrom dplyr left_join select rename '%>%' if_else distinct
#' @importFrom stats lm predict quantile sd time
#' @importFrom rlang .data
#' @importFrom igraph '%u%'
#' @export
predict_graph <- function(graphlist,
                          formulation = 2,
                          conf_level1 = NULL,
                          conf_level2 = 90,
                          dense_opt = 2,
                          weights_opt = 6,
                          weights_param = 0.001,
                          h = 1){
  # graphlist is the list of graphs
  # conf_level1 ==> this is to get different number of new nodes.
  #                 Right now, we're getting the mean predicted graph
  # conf_level2 ==> this is the upper bound of the vertex degree predictions
  #                 Right now, we're working only with the upper bound degree.
  # dense_opt   ==> this dictates of dense option needs to be used.
  #                 dense_opt = 1 for standard constraint matrix
  #                 dense_opt = 2 for dense matrix
  # weights_opt ==> this indicates the weights option
  #                 weights_opt = 1 is uniform weights
  #                 weights_opt = 2 is binary weights
  #                        1 for all except no-edges in past
  #                 weights_opt = 3 is binary weights (more selective version)
  #                        1 for all except no-edges and some 1s for new
  #                 weights_opt = 4 proportional weights

  grout_upper <- grout_lower <- graph_mean <- graph_lower <- graph_upper <- NULL

  # Step 1 - forecast the number of nodes
  pkg_message(c("i"="Predicting number of nodes"))
  new_nodes_list <-  predict_num_nodes(graphlist, conf_level1, h)
  new_nodes <- new_nodes_list$new_nodes
  if(!is.null(conf_level1)){
    # lower confidence level
    new_nodes_lower <- new_nodes_list$lower_conf
    # upper confidence level
    new_nodes_upper <- new_nodes_list$upper_conf
  }

  # Step 2 - Forecast the degree of the old nodes
  pkg_message(c("i"="Predicting old nodes degrees"))
  probj <- predict_old_nodes_degree(graphlist, conf_level2, h)

  # Step 3 - Using the above predict the mean graph
  pkg_message(c("i"="Starting internal function"))
  grout <- predict_graph_internal(graphlist,
                                  formulation,
                                  conf_level1,
                                  conf_level2,
                                  dense_opt,
                                  weights_opt,
                                  weights_param,
                                  h,
                                  probj,
                                  new_nodes)

  # Step 4 - If conf_level1 is not NULL then predict the upper and lower graphs
  if(!is.null(conf_level1)){
    # lower confidence level
    grout_lower <- predict_graph_internal(graphlist,
                                          formulation,
                                          conf_level1,
                                          conf_level2,
                                          dense_opt,
                                          weights_opt,
                                          weights_param,
                                          h,
                                          probj,
                                          new_nodes_lower)
    # upper confidence level
    grout_upper <- predict_graph_internal(graphlist,
                                          conf_level1,
                                          formulation,
                                          conf_level2,
                                          dense_opt,
                                          weights_opt,
                                          weights_param,
                                          h,
                                          probj,
                                          new_nodes_upper)
  }

  list(
    graph_mean = grout,
    graph_lower = grout_lower,
    graph_upper = grout_upper
  )
}
