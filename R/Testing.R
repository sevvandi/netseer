library(future)

setup_graph <- function(start_nodes = 1000, del_edges = 100, new_nodes = 50, new_edges = 200, num_iters = 15, mode="linear"){
  graphlist <- list()
  gr <- graphlist[[1]] <- igraph::sample_pa(start_nodes, directed = FALSE)
  if (mode == "linear"){
    graph_func <- netseer::generate_graph_linear
  } else if (mode == "exponential"){
    graph_func <- netseer::generate_graph_exp
  } else {
    stop(sprintf("UnExepected graph creation mode \"%s\", was expecting \"linear\" or \"exponential\" ", mode) )
  }
  for(i in 2:num_iters){
    gr <-  graph_func(gr, #
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
                         mode="linear", do_prof=TRUE, interval=0.01, suffix=""){
  #profile the code with Rprof, but don't try to visualise it
  #future::plan("multisession")
  fname <- paste("prof_s", start_nodes, "_nn", new_nodes, "_ae", add_edges, "_de", del_edges, suffix, sep="")
  print(fname)
  fname <- paste(fname, "Rprof", sep=".")
  print(fname)
  fname <- paste(getwd(), "profiles", fname, sep="/")
  gl <- setup_graph(start_nodes, del_edges, new_nodes, add_edges, num_iters, mode)
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

main <- function(args){
  print("Running local main")
  parser <- optparse::OptionParser()
  parser <- optparse::add_option(parser, "--start_nodes", help="Starting number of nodes", default=40)
  parser <- optparse::add_option(parser, "--add_nodes", help="Number of nodes to add for each step in the graph list", default=1)
  parser <- optparse::add_option(parser, "--add_edges", help="number of edges to add to the graph each step", default=40)
  parser <- optparse::add_option(parser, "--del_edges", help="number of edges to delete from the graph each step", default=20)
  parser <- optparse::add_option(parser, "--num_iters", help="Number of graphs to add to preceding graph list", default=15)
  parser <- optparse::add_option(parser, "--num_predict", help="Number of steps into the future to use for the prediction", default=5)
  parser <- optparse::add_option(parser, "--grow_mode", help="Whether to grow the input graphs linearly or exponentially ['linear' or 'exponential']", default="linear")
  parser <- optparse::add_option(parser, "--do_prof", action="store_true", help="Set this option to profile the code")
  parser <- optparse::add_option(parser, "--interval", help="Profiling interval", default=0.05)
  parser <- optparse::add_option(parser, "--suffix", help="suffix to add to profiling data file", default="")
  opts <- optparse::parse_args(parser, args)


  #print(commandArgs(FALSE))
  #print(opts)
  devtools::load_all()
  #(start_nodes=1000, new_nodes=30, add_edges=200, del_edges=100, num_iters=15, pred_steps=5, do_prof=TRUE, interval=0.01, suffix="", mode="linear")
  running_time <- profile_code(start_nodes=opts$start_nodes, new_nodes=opts$add_nodes, add_edges=opts$add_edges,
                               del_edges=opts$del_edges, num_iters=opts$num_iters, pred_steps=opts$num_predict,
                               mode=opts$grow_mode, do_prof=FALSE, interval=opts$interval, suffix=opts$suffix)
  print(running_time)
}

if (sys.nframe() == 0){
  main(commandArgs(TRUE))
}
