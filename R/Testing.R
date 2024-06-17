setup_graph <- function(start_nodes = 1000, del_edges = 100, new_nodes = 50, new_edges = 200, num_iters = 15){
  graphlist <- list()
  gr <- graphlist[[1]] <- igraph::sample_pa(start_nodes, directed = FALSE)
  for(i in 2:num_iters){
    gr <-  netseer::generate_graph_linear(gr, #
                              del_edge = del_edges,
                              new_nodes = new_nodes,
                              edge_increase = new_edges )
    graphlist[[i]] <- gr
  }
  return( graphlist )
}

pred_graphlist <- function(graphlist, num_steps = 5){
  predicted_graph <- netseer::predict_graph(graphlist, h=num_steps)
  return(predicted_graph)
}

profile_code <- function(start_nodes=1000, new_nodes=30, add_edges=200, del_edges=100, num_iters=15, pred_steps=5,
                         do_prof=TRUE, interval=0.01, suffix=""){
  #profile the code with Rprof, but don't try to visualise it
  fname <- paste("prof_s", start_nodes, "_nn", new_nodes, "_ae", add_edges, "_de", del_edges, suffix, sep="")
  print(fname)
  fname <- paste(fname, "Rprof", sep=".")
  print(fname)
  fname <- paste(getwd(), "profiles", fname, sep="/")
  gl <- setup_graph(start_nodes, del_edges, new_nodes, add_edges, num_iters)
  print(fname)
  start_time <- Sys.time()
  if (do_prof){
    Rprof(fname, interval = interval, line.profiling = TRUE,gc.profiling = TRUE, memory.profiling = TRUE)
  } else {
    print("Skipping profiling")
  }
  g_pred <- netseer::predict_graph(gl, h=pred_steps)
  if (do_prof){
    Rprof(NULL)
  }
  end_time <-Sys.time()
  if (do_prof){
  print(end_time - start_time)
  return(fname)
  } else {
    return (end_time - start_time)
  }
  #p <- profvis::profvis(prof_input=fname)
  #return(p)
}

#testing for the array index conversion function
quadf <- function(idx, n){
  #i <- (sqrt(1+8*idx)-1)/2
  #i <-  (1 + sqrt(-7 - 4*n + 4*n^2))/2
  #i <-  (1 + sqrt(+1 - 12*idx + 4*idx^2))/2
  #i <- (sqrt(8*idx + 4*n^2 - 12*n + 1) - 1)/2
  i <- (-sqrt(-8*idx + 4*n^2 - 4*n + 9) + 2*n + 1)/2
  #i <- as.integer(i)
  j <- n - i*(i+1)/2
  #cbind.data.frame(i, j)
  i
}

quadv <- function(n){
  vals <- 1:(n-1)
  valsum <- sum(vals)
  #quadvals <-quadf(1:valsum, valsum)
  quadvals <-quadf(1:n-1, n)
  quadvals
}

trif <- function(n){
  for(i in 1:(n-1)){
    v1 <- rep(i,(n-i))
    v2 <- (i+1):n
    df <- cbind.data.frame(v1, v2)
    if(i == 1){
      dfall <- df
    }else{
      dfall <- rbind.data.frame(dfall, df)
    }
  }
  dfall
}

#simple convert
sconv <- function(bool_string){
  idxs <- which(bool_string)
  offset <- 0
  out_idx <- 1
  n <- (sqrt(8*length(bool_string) +1) +1)/2
  src <- 1
  dst <- 1
  this_src <- n - src
  out <-matrix(data=0, nrow=length(idxs), ncol=2)
  for (idx_1d in idxs){
    while (idx_1d - offset > this_src){
      src <- src + 1
      offset <- offset + this_src
      this_src <- n - src
    }
    dst <- (idx_1d - offset) + src
    out[out_idx,1] = src
    out[out_idx,2] = dst
    out_idx <- out_idx + 1
  }
  out
}
