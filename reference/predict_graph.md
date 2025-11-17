# Predicts a graph from a time series of graphs.

This function predicts the graph at a future time step using a time
series of graphs.

## Usage

``` r
predict_graph(
  graphlist,
  formulation = 2,
  conf_level1 = NULL,
  conf_level2 = 90,
  dense_opt = 2,
  weights_opt = 8,
  weights_param = 0.001,
  h = 1
)
```

## Arguments

- graphlist:

  A list of graphs in igraph format.

- formulation:

  Formulation 2 includes an additional condition constraining total
  edges by the predicted value. Formulation 1 does not have that
  constraint. Formulation 2 gives more realistic graphs due to that
  constraint. Default is set to `2`.

- conf_level1:

  A value between 50 and 100 denoting the confidence interval for the
  number of predicted nodes in the graph. If set to `NULL` the predicted
  graph has the mean number of predicted nodes. If set to `80` for
  example, there would be 3 predicted graphs. One with mean number of
  predicted nodes, and the other two with the number of nodes
  corresponding to lower and upper confidence bounds.

- conf_level2:

  The upper confidence bound for the degree distribution. Default set to
  `90`.

- dense_opt:

  If set to `2` the dense option in R package `lpSolve` will be used.

- weights_opt:

  Weights option ranging from 1 to 6 used for different edge weight
  schemes. Weights option 1 uses uniform weights for all edges. Option 2
  uses binary weights. If the edge existed in a past graph, then weight
  is set to 1. Else set to 0. All possible new edges are assigned
  weight 1. Option 3 is a more selective version. Option 4 uses
  proportional weights according to the history. Option 5 uses
  proportional weights, but as the network is more in the past, it gives
  less weight. Option 5 uses linearly decaying proportional weights.
  Option 6 uses harmonically decaying weights. That is the network at
  `T` is given weight 1, `T-1` is given weight 1/2 and so on. Option 7
  uses 1 for edges that are present in the last graph. Option 8 is a
  slightly different to Option 7. It uses 1 for edges in the last seen
  graph and the `weights_param` for new edges. Default is set to `8`.

- weights_param:

  The weight given for possible edges from new vertices. Default set to
  `0.001`.

- h:

  The prediction time step. Default is ` h = 1`.

## Value

A list of predicted graphs. If `conf_level1` is not `NULL`, then 3
graphs are returned one with the mean number of predicted nodes and the
other 2 with the number of nodes equal to the lower and upper bound
values of prediction. If If `conf_level1` is `NULL`, only the mean
predicted graph is returned.

## Examples

``` r
set.seed(2024)
edge_increase_val <- new_nodes_val <- del_edge_val <- 0.1
graphlist <- list()
graphlist[[1]] <- gr <-  igraph::sample_pa(5, directed = FALSE)
for(i in 2:15){
  gr <-  generate_graph_exp(gr,
                        del_edge = del_edge_val,
                        new_nodes = new_nodes_val,
                        edge_increase = edge_increase_val )
  graphlist[[i]] <- gr
}
grpred <- predict_graph(graphlist[1:15], conf_level2 = 90, weights_opt = 6)
#> Registered S3 method overwritten by 'tsibble':
#>   method               from 
#>   as_tibble.grouped_df dplyr
#> Warning: 2 errors (1 unique) encountered for arima
#> [2] missing value where TRUE/FALSE needed
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo 
grpred
#> $graph_mean
#> IGRAPH 1260308 U--- 34 23 -- 
#> + edges from 1260308:
#>  [1] 28--29 17--19 11--20 11--19 10--11  9--20  5--12  5--10  3--29  3--17
#> [11]  3-- 9  3-- 6  3-- 5  2--21  2--17  2-- 4  2-- 3  1--25  1--14  1-- 6
#> [21]  1-- 5  1-- 4  1-- 2
#> 
#> $graph_lower
#> NULL
#> 
#> $graph_upper
#> NULL
#> 
```
