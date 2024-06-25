triangle_density <- function(gr){
  sum(count_triangles(gr))/(vcount(gr)*(vcount(gr)-1)*(vcount(gr)-2)/6)
}
