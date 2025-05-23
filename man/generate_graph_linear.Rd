% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_generation.R
\name{generate_graph_linear}
\alias{generate_graph_linear}
\title{Generates a bigger graph by linear growth.}
\usage{
generate_graph_linear(
  gr = NULL,
  del_edge = 1,
  new_nodes = 1,
  edge_increase = 1,
  edges_per_new_node = 3
)
}
\arguments{
\item{gr}{The input graph to generate the next graph. If set to \code{NULL}
a graph using \code{igraph::sample_pa} is used as the input graph.}

\item{del_edge}{The number of edges deleted from the input graph. Default
set to \code{1}.}

\item{new_nodes}{The number of nodes added to the input graph. Default
set to \code{1}.}

\item{edge_increase}{The number of edges added to the input graph. Default
set to \code{1}.}

\item{edges_per_new_node}{The number of edges added to the new nodes. Default
set to \code{3}.}
}
\value{
A graph.
}
\description{
Generates a bigger graph using parameters for node and edge growth. If a sequence
of graphs are created, the number of nodes would linearly increase.
}
\examples{
set.seed(1)
gr <- generate_graph_linear()
gr

}
