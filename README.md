
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Predicting graph structure from a time series of graphs <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/sevvandi/netseer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sevvandi/netseer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

If you have a time series of dynamic graphs with changing structure, how
would you predict future graphs? This is the goal of netseer.

Netseer predicts the graph structure including new nodes and edges from
a time series of graphs. It adapts Flux Balance Analysis, a method used
in metabolic network reconstruction to predict the structue of future
graphs. The methodology is explained in the preprint (Kandanaarachchi et
al. 2025).

## R package

The algorithm is available in both R and Python. The R package `netseer`
in on CRAN and can be installed as follows:

``` r
install_packages("netseer")
```

The vignette for the R package is available under [Get
Started](articles/netseer.html) at <https://sevvandi.github.io/netseer/>

## Python code

The Python code is available at
<https://github.com/sevvandi/netseer-python/>. The Python package will
be made available in the near future.

## Coding Credits

A big shout out to Stefan Westerlund and Brodie Oldfield for helping
with this package. Stefan optimized the algorithm in R and coded it from
scratch in Python. Brodie is currently packaging up the Python code and
is working on a website for the Python version.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-kand2025graphpred" class="csl-entry">

Kandanaarachchi, Sevvandi, Ziqi Xu, Stefan Westerlund, and Conrad
Sanderson. 2025. “Predicting Graph Structure via Adapted Flux Balance
Analysis.” <https://arxiv.org/abs/2507.05806>.

</div>

</div>
