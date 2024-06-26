#' Generates a bigger graph using exponential growth.
#'
#' Generates a bigger graph using parameters for node and edge growth. If a sequence
#' of graphs are created, the number of nodes in this sequence would exponentially increase.
#'
#'@param gr The input graph to generate the next graph. If set to \code{NULL}
#'a graph using \code{igraph::sample_pa} is used as the input graph.
#'@param del_edge The proportion of edges deleted from the input graph. Default
#'set to \code{0.1}.
#'@param new_nodes The proportion of nodes added to the input graph. Default
#'set to \code{0.1}.
#'@param edge_increase The proportion of edges added to the input graph. Default
#'set to \code{0.1}.
#'
#'@return A graph.
#'
#'@examples
#'set.seed(1)
#'gr <- generate_graph_exp()
#'gr
#'
#'@export
generate_graph_exp <- function(gr = NULL, del_edge = 0.1, new_nodes = 0.1, edge_increase = 0.1){
  # gr - graph to start with
  # del_edge - between 0 and 1. The proportion of edges to delete
  # new_nodes - if less than 1 then it is a proportion, else it is the number of nodes to add
  # edge_increase - the proportion of edges to add

  if(is.null(gr)){
    gr <- igraph::sample_pa(5, directed = FALSE)
  }

  if(new_nodes < 1){
    new_nodes <- ceiling(igraph::vcount(gr)*(new_nodes))
  }
  num_edges <- igraph::ecount(gr)
  edges <- igraph::E(gr)

  # Removing edges
  num_edges1 <- floor(del_edge*num_edges)
  edges_remove <- sample(igraph::E(gr), num_edges1)
  gr2 <- igraph::delete_edges(gr, edges_remove)

  # Adding vertices
  gr2 <- igraph::add_vertices(gr2, new_nodes)

  # Add edges
  num_edges2 <- ceiling(num_edges*(1 + edge_increase)) -  igraph::ecount(gr2)
  gr2 <- gr2 + igraph::edge(sample(igraph::V(gr2), num_edges2*2, replace = T))
  gr2 <- igraph::simplify(gr2)
  return(gr2)
}


#' Generates a bigger graph by linear growth.
#'
#' Generates a bigger graph using parameters for node and edge growth. If a sequence
#' of graphs are created, the number of nodes would linearly increase.
#'
#'@param gr The input graph to generate the next graph. If set to \code{NULL}
#'a graph using \code{igraph::sample_pa} is used as the input graph.
#'@param del_edge The number of edges deleted from the input graph. Default
#'set to \code{1}.
#'@param new_nodes The number of nodes added to the input graph. Default
#'set to \code{1}.
#'@param edge_increase The number of edges added to the input graph. Default
#'set to \code{1}.
#'
#'@return A graph.
#'
#'@examples
#'set.seed(1)
#'gr <- generate_graph_linear()
#'gr
#'
#'@export
#'@export
generate_graph_linear <- function(gr = NULL, del_edge = 1, new_nodes = 1, edge_increase = 1){

  if(is.null(gr)){
    gr <- igraph::sample_pa(10, directed = FALSE)
  }


  num_edges <- igraph::ecount(gr)
  edges <- igraph::E(gr)

  # Removing edges
  edges_remove <- sample(igraph::E(gr), del_edge)
  gr <- igraph::delete_edges(gr, edges_remove)

  # Possible edges to add
  adj1 <- igraph::as_adjacency_matrix(gr)
  adj2 <- adj1 %*%adj1
  possible_edges <- c()
  for(ll in 1:(NROW(adj2)-1)){
    for(mm in (ll+1):NCOL(adj2)){
      if((adj2[ll, mm] > 0) & (adj1[ll, mm] == 0)){
        possible_edges <- c(possible_edges, c(ll, mm))
      }
    }
  }

  # Add edges
  num_edges2 <- edge_increase*2
  edges_to_add_inds <- sort(sample(length(possible_edges)/2, num_edges2))
  edges_to_add_inds1 <- edges_to_add_inds*2 -1
  edges_to_add_inds2 <- edges_to_add_inds*2
  edges_to_add <- possible_edges[c(rbind(edges_to_add_inds1, edges_to_add_inds2))]
  gr2 <- gr + igraph::edge(edges_to_add)

  # Adding vertices
  gr2 <- igraph::simplify(gr2)
  gr3 <- igraph::add_vertices(gr2, new_nodes)

  # Add edges to these vertices
  num_edges3 <- edge_increase*3
  probs <- igraph::degree(gr2)/sum(igraph::degree(gr2))
  new_node_ids <-  which(igraph::degree(gr3) == 0)
  e1 <- sample(igraph::V(gr2), num_edges3, replace = T, prob = probs )
  e2 <- sample(new_node_ids, num_edges3, replace = T)
  gr3 <- gr3 + igraph::edge(rbind(e1, e2))
  gr3 <- igraph::simplify(gr3)


  return(gr3)
}
