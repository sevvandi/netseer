# Reads graphs from a folder

This function reads graphs from a folder to a list

## Usage

``` r
read_graph_list(path_to_graphs, format)
```

## Arguments

- path_to_graphs:

  The folder where graphs are contained.

- format:

  Formats supported by
  [`igraph::read_graph`](https://r.igraph.org/reference/read_graph.html).

## Value

A list of graphs in `igraph` format.

## Examples

``` r
# example code
path_to_graphs <- system.file("extdata", package = "netseer")
grlist <- read_graph_list(path_to_graphs = path_to_graphs, format = "gml")
grlist
#> [[1]]
#> IGRAPH b959539 U--- 46 116 -- Barabasi graph
#> + attr: name (g/c), power (g/n), m (g/n), zeroappeal (g/n), algorithm
#> | (g/c), id (v/n)
#> + edges from b959539:
#>  [1]  1-- 2  1-- 4  1-- 5  1-- 9  1--11  1--12  1--14  1--16  1--17  1--20
#> [11]  1--22  1--33  1--40  1--43  2-- 4  2--10  2--11  2--13  2--14  2--15
#> [21]  2--19  2--22  2--28  2--29  2--38  2--39  3-- 5  3-- 6  3--11  3--13
#> [31]  3--14  3--21  3--28  3--33  3--35  3--37  3--41  3--46  4-- 8  4-- 9
#> [41]  4--10  4--14  4--16  4--17  4--18  4--21  4--22  4--30  4--32  4--36
#> [51]  4--43  5-- 6  5-- 7  5-- 9  5--11  5--17  5--18  5--23  5--25  5--30
#> [61]  5--44  6--18  6--19  6--20  7--25  7--34  7--43  9--31  9--37 10--13
#> + ... omitted several edges
#> 
#> [[2]]
#> IGRAPH ee98131 U--- 51 130 -- Barabasi graph
#> + attr: name (g/c), power (g/n), m (g/n), zeroappeal (g/n), algorithm
#> | (g/c), id (v/n)
#> + edges from ee98131:
#>  [1]  1-- 2  1-- 4  1-- 5  1-- 9  1--11  1--12  1--14  1--16  1--17  1--20
#> [11]  1--22  1--33  1--40  1--43  2-- 4  2--10  2--11  2--13  2--14  2--15
#> [21]  2--19  2--22  2--28  2--29  2--38  2--39  3-- 5  3-- 6  3--11  3--13
#> [31]  3--14  3--21  3--28  3--33  3--35  3--37  3--41  3--46  3--48  4-- 8
#> [41]  4-- 9  4--10  4--14  4--16  4--17  4--18  4--21  4--22  4--30  4--32
#> [51]  4--36  4--43  5-- 6  5-- 7  5-- 9  5--11  5--17  5--18  5--23  5--25
#> [61]  5--30  5--44  5--50  6--18  6--19  6--20  6--35  6--47  7--25  7--34
#> + ... omitted several edges
#> 
#> [[3]]
#> IGRAPH 9322b56 U--- 56 145 -- Barabasi graph
#> + attr: name (g/c), power (g/n), m (g/n), zeroappeal (g/n), algorithm
#> | (g/c), id (v/n)
#> + edges from 9322b56:
#>  [1] 1-- 2 1-- 4 1-- 5 1-- 9 1--11 1--12 1--14 1--16 1--17 1--20 1--22 1--32
#> [13] 1--33 1--40 1--43 2-- 4 2--10 2--11 2--13 2--14 2--15 2--19 2--22 2--28
#> [25] 2--29 2--38 2--39 2--56 3-- 5 3-- 6 3--11 3--13 3--14 3--21 3--28 3--33
#> [37] 3--35 3--37 3--41 3--46 3--48 4-- 8 4-- 9 4--10 4--16 4--17 4--18 4--21
#> [49] 4--22 4--30 4--32 4--36 4--43 4--53 5-- 6 5-- 7 5-- 9 5--11 5--17 5--18
#> [61] 5--23 5--25 5--30 5--44 5--50 5--56 6--18 6--19 6--20 6--35 6--47 7--25
#> [73] 7--34 7--43 7--52 9--31 9--37 9--56
#> + ... omitted several edges
#> 
```
