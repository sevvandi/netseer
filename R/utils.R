triangle_density <- function(gr){
  sum(igraph::count_triangles(gr))/(igraph::vcount(gr)*(igraph::vcount(gr)-1)*(igraph::vcount(gr)-2)/6)
}
