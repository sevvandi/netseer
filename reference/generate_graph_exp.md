# Generates a bigger graph using exponential growth.

Generates a bigger graph using parameters for node and edge growth. If a
sequence of graphs are created, the number of nodes in this sequence
would exponentially increase.

## Usage

``` r
generate_graph_exp(
  gr = NULL,
  del_edge = 0.1,
  new_nodes = 0.1,
  edge_increase = 0.1
)
```

## Arguments

- gr:

  The input graph to generate the next graph. If set to `NULL` a graph
  using
  [`igraph::sample_pa`](https://r.igraph.org/reference/sample_pa.html)
  is used as the input graph.

- del_edge:

  The proportion of edges deleted from the input graph. Default set to
  `0.1`.

- new_nodes:

  The proportion of nodes added to the input graph. Default set to
  `0.1`.

- edge_increase:

  The proportion of edges added to the input graph. Default set to
  `0.1`.

## Value

A graph.

## Examples

``` r
set.seed(1)
gr <- generate_graph_exp()
gr
#> IGRAPH 76687d7 U--- 6 4 -- Barabasi graph
#> + attr: name (g/c), power (g/n), m (g/n), zero.appeal (g/n), algorithm
#> | (g/c)
#> + edges from 76687d7:
#> [1] 1--2 1--3 1--4 2--5
```
