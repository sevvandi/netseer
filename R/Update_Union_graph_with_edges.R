# --------------------------------------------------------------------
# EXAMPLE 1 : Update the graph with the edges that are not present in the original graph
# EXAMPLE 2 : Evaluate predicted graph
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# EXAMPLE 1 : Update the graph with the edges that are not present in the original graph
# --------------------------------------------------------------------
#library(igraph)
add_edges_exampple <- function(){
  set.seed(1)
  biggr1 <- biggr <- erdos.renyi.game(10, 0.2, directed = FALSE)
  plot(biggr)
  A <- igraph::as_adjacency_matrix(biggr)
  Asq <- A %*% A
  for(ii in 1:(NROW(Asq) - 1)){
    for(jj in (ii + 1):NROW(Asq)){
      if(Asq[ii, jj] > 0){
        if(A[ii, jj] == 0){
          biggr <- igraph::add_edges(biggr, c(ii, jj), weight = 0.05)
        }
      }
    }
  }

  par(mfrow = c(1,2))
  plot(biggr1)
  E(biggr)$weight[E(biggr)$weight == 0] <- 1
  plot(biggr, edge.label = E(biggr)$weight)
}


# --------------------------------------------------------------------
# EXAMPLE 2 : Evaluate predicted graph
# --------------------------------------------------------------------

# Function2 to evaluate graph prediction

eval_metrics <- function(actualgr, predgr){

  new_vert <- NULL
  metric1 <- (vcount(predgr) - vcount(actualgr))/vcount(actualgr)

  if(metric1 > 0){
    # vcount(predgr) > vcount(actualgr)
    new_vert <- vcount(predgr) - vcount(actualgr)
    actualgr <-igraph::add_vertices(actualgr,nv =new_vert)
  }else if(metric1 < 0){
    # vcount(predgr) < vcount(actualgr)
    new_vert <-  vcount(actualgr) - vcount(predgr)
    predgr <-igraph::add_vertices(predgr,nv =new_vert)
  }
  adj_act <-igraph::as_adj(actualgr)
  adj_pred <- igraph::as_adj(predgr)
  metric2 <- (Matrix::norm((adj_pred - adj_act), type = "F"))/(Matrix::norm(adj_act, type = "F"))

  metric3 <- diff_metrics(as.vector(adj_act), as.vector(adj_pred))

  deg_act <- igraph::degree(actualgr)
  deg_prd <- igraph::degree(predgr)
  deg_cosine_sim <- (sum(deg_act*deg_prd))/(norm(as.matrix(deg_act), type = "F")*norm(as.matrix(deg_prd), type = "F"))
  # lsa::cosine(deg_act, deg_prd)

  lap_act <- igraph::laplacian_matrix(actualgr)
  lap_prd <- igraph::laplacian_matrix(predgr)

  eig_act <- eigen(lap_act)$values
  eig_prd <- eigen(lap_prd)$values
  eig_cosine_sim <- (sum(eig_act*eig_prd))/(norm(as.matrix(eig_act), type = "F")*norm(as.matrix(eig_prd), type = "F"))


  list(
    node_error = metric1,
    adjacency_error = metric2,
    degree_cosine = deg_cosine_sim,
    lap_eigen_cosine = eig_cosine_sim,
    precision = metric3$precision,
    recall = metric3$recall,
    fmeasure = metric3$fmeasure,
    sensitivity = metric3$sensitivity,
    specificity = metric3$specificity,
    gmean = metric3$gmean,
    accuracy = metric3$accuracy

  )
}


# Evaluate prediction function
diff_metrics <- function(act, pred) {
  # positives to be denoted by 1 and negatives with 0
  stopifnot(length(act) == length(pred))
  n <- length(act)
  tp <- sum((act == 1) & (pred == 1))
  tn <- sum((act == 0) & (pred == 0))
  fp <- sum((act == 0) & (pred == 1))
  fn <- sum((act == 1) & (pred == 0))
  prec <- (tp + tn) / n
  sn <- tp / (tp + fn)
  sp <- tn / (tn + fp)
  precision <- if_else(
    (tp + fp) == 0,
    0,
    tp / (tp + fp)
  )
  recall <- tp / (tp + fn)
  fmeasure <- if_else(
    (precision == 0) & (recall == 0),
    0,
    2 * precision * recall / (precision + recall)
  )
  tibble::tibble(
    N = n,
    true_pos = tp,
    true_neg = tn,
    false_pos = fp,
    false_neg = fn,
    accuracy = prec,
    sensitivity = sn,
    specificity = sp,
    gmean = sqrt(sn * sp),
    precision = precision,
    recall = recall,
    fmeasure = fmeasure
  )
}

example_test <-function(){
  # Example given in the function
  #library(igraph)
  #library(netseer)
  #library(tidyverse)
  set.seed(2024)
  edge_increase_val <- new_nodes_val <- del_edge_val <- 0.1
  graphlist <- list()
  graphlist[[1]] <- gr <-  igraph::sample_pa(5, directed = FALSE)
  for(i in 2:16){
    gr <-  netseer::generate_graph_exp(gr,
                              del_edge = del_edge_val,
                              new_nodes = new_nodes_val,
                              edge_increase = edge_increase_val )
    graphlist[[i]] <- gr
  }
  grpred <- predict_graph(graphlist[1:15], conf_level2 = 90, weights_opt = 6)

  eval_metrics(graphlist[[16]], grpred$graph_mean)
}




