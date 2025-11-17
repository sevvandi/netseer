# Gives an error measurement for predicted graphs

This function compares the predicted graph with the actual and comptues
the node and edge error as a proportion

## Usage

``` r
measure_error(actual, predicted)
```

## Arguments

- actual:

  The ground truth or actual graph.

- predicted:

  The predicted graph.

## Value

The node error and edge error as a proportion.

## Examples

``` r
data(syngraphs)
# Taking the 20th graph as the actual and the 19th graph as predicted.
measure_error(syngraphs[[20]], syngraphs[[19]])
#> $node_error
#> [1] 0.08928571
#> 
#> $edge_error
#> [1] 0.1034483
#> 

```
