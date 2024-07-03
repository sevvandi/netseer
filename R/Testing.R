library(future)

setup_graph <- function(start_nodes = 1000, del_edges = 100, new_nodes = 50, new_edges = 200, nn_edges=3, num_iters = 15, mode="linear"){
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
                              edge_increase = new_edges,
                              edges_per_new_node = nn_edges
                      )
    graphlist[[i]] <- gr
  }
  return( graphlist )
}

pred_graphlist <- function(graphlist, num_steps = 5){
  predicted_graph <- netseer::predict_graph(graphlist, h=num_steps)
  return(predicted_graph)
}

profile_code <- function(start_nodes=1000, new_nodes=30, add_edges=200, del_edges=100, nn_edges=3, num_iters=15, pred_steps=5,
                         mode="linear", do_prof=TRUE, interval=0.05, suffix="", verbosity="verbose"){
  #profile the code with Rprof, but don't try to visualise it
  options(netseer.verbose=verbosity)
  pkg_message(c("i"="Generating random graph list"))
  gl <- setup_graph(start_nodes, del_edges, new_nodes, add_edges, nn_edges, num_iters, mode)
  start_time <- Sys.time()
  if (do_prof){
    fname <- paste("prof_s", start_nodes, "_nn", new_nodes, "_ae", add_edges, "_de", del_edges, suffix, sep="")
    fname <- paste(fname, "Rprof", sep=".")
    fname <- paste(getwd(), "profiles", fname, sep="/")
    pkg_message(c("i"=sprint("Saving profiling results to %s",fname)))
    Rprof(fname, interval = interval, line.profiling = TRUE,gc.profiling = TRUE, memory.profiling = TRUE)
  } else {
    pkg_message(c("i"="Skipping profiling"))
  }
  pkg_message(c("i"="Running graph prediction"))
  g_pred <- netseer::predict_graph(gl, h=pred_steps)
  if (do_prof){
    Rprof(NULL)
  }
  end_time <-Sys.time()
  if (do_prof){
    pkg_message(toString(end_time - start_time))
  return(fname)
  } else {
    return (end_time - start_time)
  }
}

main <- function(args){
  parser <- optparse::OptionParser()
  parser <- optparse::add_option(parser, "--start_nodes", help="Starting number of nodes", default=40)
  parser <- optparse::add_option(parser, "--add_nodes", help="Number of nodes to add for each step in the graph list", default=1)
  parser <- optparse::add_option(parser, "--add_edges", help="number of edges to add to the graph each step", default=40)
  parser <- optparse::add_option(parser, "--del_edges", help="number of edges to delete from the graph each step", default=20)
  parser <- optparse::add_option(parser, "--nn_edges", help="number of edges to add per newly-added node", default=3)
  parser <- optparse::add_option(parser, "--num_iters", help="Number of graphs to add to preceding graph list", default=15)
  parser <- optparse::add_option(parser, "--num_predict", help="Number of steps into the future to use for the prediction", default=5)
  parser <- optparse::add_option(parser, "--grow_mode", help="Whether to grow the input graphs linearly or exponentially ['linear' or 'exponential']", default="linear")
  parser <- optparse::add_option(parser, "--do_prof", action="store_true", help="Set this option to profile the code")
  parser <- optparse::add_option(parser, "--interval", help="Profiling interval", default=0.05)
  parser <- optparse::add_option(parser, "--suffix", help="suffix to add to profiling data file", default="")
  parser <- optparse::add_option(parser, "--verbosity", help="Verbosity level - how much information to print ['quiet', 'verbose' or 'debug']", default="quiet")
  opts <- optparse::parse_args(parser, args)


  #print(commandArgs(FALSE))
  #print(opts)
  devtools::load_all()
  #(start_nodes=1000, new_nodes=30, add_edges=200, del_edges=100, num_iters=15, pred_steps=5, do_prof=TRUE, interval=0.01, suffix="", mode="linear")
  running_time <- profile_code(start_nodes=opts$start_nodes, new_nodes=opts$add_nodes, add_edges=opts$add_edges,
                               del_edges=opts$del_edges, nn_edges=opts$nn_edges, num_iters=opts$num_iters, pred_steps=opts$num_predict,
                               mode=opts$grow_mode, do_prof=FALSE, interval=opts$interval, suffix=opts$suffix, verbosity=opts$verbosity)
  print(running_time)
}

if (sys.nframe() == 0){
  main(commandArgs(TRUE))
}
