library(future)
# THIS FILE DOES THE FOLLOWING
# OPTIMIZATION PROBLEM = 1st ONE (0 \leq Su \leq f_u(d)
# IMPLEMENTATION DETAILS: dense matrix in lpSolve and
#                         standard solving in lpSolve
# WEIGHT SCHEMES:given in weights_opt
# PREDICTION/FORECASTING STEP: h in T+h
# PREDICTED NUMBER OF NODES: conf_level1
#                            (only consider increasing in size graphs)
# CONSTRAINTS: PREVIOUS - one constraint for degree of each node
#                       - one additional constraint for total number of edges
# OPTIMIZATION: Make the optimization tighter.
#               Remove columns with zero weights in the big matrix
#               because the zero weights mean no edge.
#               This will make a smaller optimization problem.
#               Remove rows with zero rhs.
#




predict_graph_internal <- function(graphlist,
                                   formulation,
                                   conf_level1 = NULL,
                                   conf_level2 = 90,
                                   weights_opt = 2,
                                   weights_param = 0.25,
                                   h = 1,
                                   probj,
                                   new_nodes){
  # graphlist is the list of graphs
  # conf_level1 ==> this is to get different number of new nodes.
  #                 Right now, we're getting the mean predicted graph
  # conf_level2 ==> this is the upper bound of the vertex degree predictions
  #                 Right now, we're working only with the upper bound degree.
  # dense_opt   ==> this dictates of dense option needs to be used.
  #                 dense_opt = 1 for standard constraint matrix
  #                 dense_opt = 2 for dense matrix
  # weights_opt ==> this indicates the weights option
  #                 weights_opt = 1 is uniform weights
  #                 weights_opt = 2 is binary weights
  #                        1 for all except no-edges in past
  #                 weights_opt = 3 is binary weights (more selective version)
  #                        1 for all except no-edges and some 1s for new
  #                 weights_opt = 4 proportional weights

  pkg_debug(c("i"=sprintf("Using weights option %d", weights_opt)))
  pkg_message(c("i"="Setting maximum degrees constraint"))
  degree_constraints <- probj$degree_hi #the maximum degrees for each vertex in the predicted graph
  total_edges_constraint <- probj$total_edges_upper
  if (sum(new_nodes) > 0){
    #determine the degree(s) that the new nodes should have
    new_nodes_degrees_per_age <- get_new_nodes_degrees(graphlist, h, conf_level2)
    new_node_degree_constraints <- new_nodes_edge_constraints(new_nodes, new_nodes_degrees_per_age)
    num_new_node_edges <- sum(new_node_degree_constraints)
    total_edges_constraint <- total_edges_constraint - num_new_node_edges
  }
  new_nodes <- sum(new_nodes) #set new nodes from nodes per predicted timestep to total number of new nodes
  pkg_message(c("i"="Constructing union graph"))
  biggr <- construct_union_graph(graphlist, new_nodes, degree_constraints, rm_verts, weights_opt = weights_opt, weights_param = weights_param)
  pkg_debug(c("i"=sprintf("The union graph has %d nodes, %d edges, with constraint to use %d edges", igraph::vcount(biggr), igraph::ecount(biggr), total_edges_constraint)))
  pkg_message(c("i"="Setting sparse solver constraints"))
  lpsolve_vars <- sparse_solver_formulation(biggr, degree_constraints, total_edges_constraint)
  pkg_message(c("i"="Running LPSolve"))
  lpobj <- run_lpsolve_dense(lpsolve_vars$constraint_matrix, lpsolve_vars$objective_function, lpsolve_vars$constraint_rhs, lpsolve_vars$constraint_direction)
  pkg_debug(c("i"=sprintf("LPSolve solution status (should be 0) is %d", lpobj$status) ))
  pkg_message(c("i"="Constructing graph from LPSove Solution"))
  grout <- graph_from_solution(lpobj$solution, lpsolve_vars$used_edges, igraph::vcount(biggr))
  if (new_nodes > 0){
    pkg_message(c("i"="Adding new nodes and their edges"))
    grout <- add_new_nodes_and_edges(grout, new_node_degree_constraints, degree_constraints)
  }


  if (getOption("netseer.verbose", "quiet") == "debug"){
    pkg_debug(c("i"="Checking that constraints are met"))
    if (new_nodes > 0){
      maximum_degrees <- c(degree_constraints, new_node_degree_constraints)
    }  else {
      maximum_degrees <- degree_constraints
    }
    check_constraints(grout, maximum_degrees, probj$total_edges_upper)
  }

  pkg_message(c("i"="Predict graph internal done"))
  grout

}

predict_num_nodes <- function(graphlist, conf_level = 90, h = 1){
  # Forecast the number of nodes
  NN <- length(graphlist)
  num_nodes <- rep(0, NN)
  nodes <- index <- lower_conf <- upper_conf <- NULL
  for(kk in 1:NN){
    num_nodes[kk] <- igraph::gorder(graphlist[[kk]])
  }

  df <- data.frame(index = 1:NN, nodes = num_nodes)

  if(NN < 15){
    lmmod <- lm(nodes ~index, data = df)
    next_nodes <- predict(lmmod, newdata = data.frame(index = (NN+h)))
  }else{
    next_nodes <- df %>% tsibble::as_tsibble(index = index) %>%
      fabletools::model(ets = fable::ETS(nodes)) %>%
      fabletools::forecast(h = h) %>%
      pull(.data$.mean)
    fit <- df %>% tsibble::as_tsibble(index = index) %>%
      fabletools::model(ets = fable::ETS(nodes)) %>%
      fabletools::forecast(h = h)
    if(!is.null(conf_level)){
      dfhilo <- fit %>%
        fabletools::hilo(level = conf_level)
      lower_conf <-  floor(unlist(dfhilo[ ,5])[1]) - igraph::gorder(graphlist[[NN]])
      upper_conf <-  ceiling(unlist(dfhilo[ ,5])[2]) - igraph::gorder(graphlist[[NN]])
    }
  }


  new_nodes <-  ceiling(next_nodes) - igraph::gorder(graphlist[[NN]])
  new_nodes <- c(new_nodes[[1]], diff(new_nodes)) #new nodes for each timestep
  #handle the possibility of the number of nodes being reduces at some time steps, and being increased at later steps
  #we prefer to handle "removed" nodes as nodes with zero edges attached
  removed_nodes <- 0
  for (i in 1:length(new_nodes)){
    if (new_nodes[[i]] < 0){
      removed_nodes <- removed_nodes + abs(new_nodes[[i]])
      new_nodes[[i]] <- 0
      next
    }
    if (removed_nodes <= 0) next
    this_removed_nodes <- min(removed_nodes, new_nodes[[i]])
    removed_nodes <- removed_nodes - this_removed_nodes
    new_nodes[[i]] <- new_nodes[[i]] - this_removed_nodes
  }
  if (removed_nodes > 0){
    new_nodes[[length(new_nodes)]] <- new_nodes[[length(new_nodes)]] - removed_nodes
  }
  list(
    new_nodes = new_nodes,#[h]
    lower_conf = lower_conf,
    upper_conf = upper_conf
  )

}


construct_union_graph <- function(graphlist,
                                  new_nodes,
                                  new_deg,
                                  rm_verts = NULL,
                                  weights_opt = 2,
                                  weights_param = 0.25){
  # Construct the union graph
  # new_nodes = number of new nodes
  # new_deg =  average degree of new nodes
  # rm_verts = if new nodes is negative, then the set of nodes to remove
  # weights_opt = 2: Binary weights - new nodes connected to all old nodes
  # weights_opt = 3: Binary weights - new nodes connected to most connected old nodes
  # weights_opt = 4: Proportional weights

  NN <- length(graphlist)

  if(weights_opt == 4){
    # Graph needs weights
    for(ii in 1:NN){
      grr <- graphlist[[ii]]
      igraph::E(grr)$weight <- 1
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
        igraph::E(biggr)$weight_1[is.na(igraph::E(biggr)$weight_1)] <- 0
        igraph::E(biggr)$weight_2[is.na(igraph::E(biggr)$weight_2)] <- 0
        igraph::E(biggr)$weight <- igraph::E(biggr)$weight_1 +  igraph::E(biggr)$weight_2
      }
    }
    igraph::E(biggr)$weight <- igraph::E(biggr)$weight/max(igraph::E(biggr)$weight)
  }else if(weights_opt == 5){
    # This is linearly decaying weights over time.
    # ii goes from 1 to NN
    # ii = 1 is the first time step
    # ii = NN is the most recent time step
    for(ii in 1:NN){
      grr <- graphlist[[ii]]
      igraph::E(grr)$weight <- ii # this line make linear decaying weights
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
        igraph::E(biggr)$weight_1[is.na(igraph::E(biggr)$weight_1)] <- 0
        igraph::E(biggr)$weight_2[is.na(igraph::E(biggr)$weight_2)] <- 0
        igraph::E(biggr)$weight <- igraph::E(biggr)$weight_1 +  igraph::E(biggr)$weight_2
      }
    }
    igraph::E(biggr)$weight <- igraph::E(biggr)$weight/max(igraph::E(biggr)$weight)
  }else if(weights_opt == 6){
    # This is decaying weights over time like 1, 1/2, 1/3, etc . . .
    # ii goes from 1 to NN
    # ii = 1 is the first time step
    # ii = NN is the most recent time step
    for(ii in 1:NN){
      grr <- graphlist[[ii]]
      igraph::E(grr)$weight <- 1/(NN - ii + 1) # this line make harmonically decaying weights
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
        igraph::E(biggr)$weight_1[is.na(igraph::E(biggr)$weight_1)] <- 0
        igraph::E(biggr)$weight_2[is.na(igraph::E(biggr)$weight_2)] <- 0
        igraph::E(biggr)$weight <- igraph::E(biggr)$weight_1 +  igraph::E(biggr)$weight_2
      }
    }
    igraph::E(biggr)$weight <- igraph::E(biggr)$weight/max(igraph::E(biggr)$weight)
  }else if(weights_opt == 7){
    # This set of weights give 1 to the last time step and 0 to the rest
    biggr <- graphlist[[NN]]
    igraph::E(biggr)$weight <- 1
  }else{
  # For weight options 1, 2, 3
  # For t = 1 to n
    for(ii in 1:NN){
      grr <- graphlist[[ii]]
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
      }
    }
  }
  existing_edges <- igraph::ecount(biggr)

  # Weights for new edges according to quantile
  if(weights_opt %in% c(4,5,6)){
    new_weights <- quantile(igraph::E(biggr)$weight, probs = weights_param)
  }else if(weights_opt == 7){
    new_weights <- weights_param
  }


  # Add edges if nodes share a neighbour but do not already have an edge directly to each other
  pkg_message(c("i"="Adding neighbours of neighbours"))
  bigadj1 <- igraph::as_adjacency_matrix(biggr)
  bigadj2 <- bigadj1 %*% bigadj1
  non_edges <- Matrix::which((bigadj2 > 0) & (bigadj1 == 0), arr.ind=TRUE) #neighbours of neighbours edges
  non_edges <- t(non_edges[non_edges[,1] < non_edges[,2],])

  if(weights_opt %in% c(4,5,6,7)){
    biggr <- igraph::add_edges(biggr, non_edges, weight = 2*new_weights)
  }else{
    # weights_opt 1, 2, 3
    biggr <- igraph::add_edges(biggr, non_edges)
  }

  # if(new_nodes > 0){
  #   # Add new nodes
  #   pkg_message(c("i"="Adding new nodes, and their edges"))
  #   biggr <- igraph::add_vertices(biggr, new_nodes)
  #   new_vertices <- igraph::V(biggr)[( igraph::gorder(biggr)-new_nodes+1):igraph::gorder(biggr)]
  #   # Use only a fixed number of old vertices
  #   # Which vertices have the highest degree
  #   grlast <- graphlist[[NN]]
  #   num_attach <- 10 # attach potential edges to this number of nodes in the union graph, with the highest degrees
  #   pkg_debug(c("i"=sprintf("Attaching possible edges to the new nodes and the top %d existing nodes, by degree", num_attach)))
  #   if(igraph::vcount(grlast) > num_attach){
  #     old_vertices <- igraph::V(biggr)[order(igraph::degree(grlast), decreasing = TRUE)[1:num_attach]]
  #   }else{
  #     old_vertices <- igraph::V(biggr)[1:(igraph::gorder(biggr)-new_nodes)]
  #   }
  #
  #   new_deg2 <- min(2*new_deg, igraph::vcount(biggr))
  #   verts <- order(igraph::degree(biggr), decreasing = TRUE)[1:new_deg2]
  #
  #   if((weights_opt == 1)|(weights_opt == 2)){
  #     # Binary weights - new nodes connected to all old nodes
  #     # Add new edges from new nodes to all old nodes
  #     possible_edges <- c(rbind(rep(old_vertices, new_nodes), rep(new_vertices, each = length(old_vertices)) ))
  #   }else if(weights_opt == 3){
  #     # Binary weights - new nodes connected to most connected old nodes
  #     # Add new edges from new nodes to mostly connected vertices
  #     possible_edges <- c(rbind(rep(verts, new_nodes), rep(new_vertices, each = length(verts)) ))
  #   }else if(weights_opt == 4|weights_opt == 5|weights_opt == 6|weights_opt == 7){
  #     # Proportional weights - new nodes connected to all old nodes
  #     # But the weights will be much smaller
  #     possible_edges <- c(rbind(rep(old_vertices, new_nodes), rep(new_vertices, each = length(old_vertices)) ))
  #     new_weights0 <- igraph::degree(biggr)[old_vertices]/(sum(igraph::degree(biggr)[old_vertices]))
  #     new_weights <- quantile(igraph::E(biggr)$weight, probs = new_weights0)
  #     new_weights <- rep(new_weights, new_nodes )
  #   }
  #   if(weights_opt == 4|weights_opt == 5|weights_opt == 6|weights_opt == 7){
  #     biggr <- biggr %>%
  #       igraph::add_edges(possible_edges,weight = new_weights)
  #   }else{
  #     biggr <- biggr %>%
  #       igraph::add_edges(possible_edges)
  #   }
  #
  # }

  total_edges <- igraph::ecount(biggr)
  new_edges <- total_edges - existing_edges
  pkg_debug(c("i"=sprintf("The union graph has %d existing edges, %d new edges, for %d total edges and %.2f%% are new nodes", existing_edges, new_edges, total_edges, 100*new_edges/total_edges)))

  biggr
}


predict_old_nodes_degree <- function(graphlist, conf_level2, h){
  # Forecast the degree of the old nodes
  future::plan(future::multisession)
  degree <- edges <- hilow <- lower <- upper <- upper2 <- lower2  <- mean2 <- NULL

  NN <- length(graphlist)
  for(jj in 1:NN){
    gr <- graphlist[[jj]]
    if(is.null(names(igraph::V(gr)))){
      v_names <- unclass(igraph::V(gr))
    }else{
      v_names <- names(igraph::V(gr))
    }
    df <- data.frame(time = jj, vertex = v_names, degree = igraph::degree(gr))
    df_sum <- data.frame(time = jj, edges = igraph::ecount(gr))
    if(jj == 1){
      dfall <- df
      dfall_sum <- df_sum
    }else{
      dfall <- rbind.data.frame(dfall, df)
      dfall_sum <- rbind.data.frame(dfall_sum, df_sum)
    }
  }

  dfmerge <- tibble::as_tibble(dfall) %>%
    tsibble::as_tsibble(index = time , key = vertex) %>%
    dplyr::arrange(time, vertex) %>%
    tsibble::fill_gaps(degree = 0)

  dffreq <- dfall %>%
    group_by(vertex) %>%
    summarize(freq = n())


  # Fit ARIMA models
  #original version
  #fit <- dfmerge %>%
  #  fabletools::model(arima = fable::ARIMA(degree),
  #                    naive = fable::NAIVE(degree)) %>%
  #  fabletools::forecast(h = h) %>%
  #  full_join(dffreq)

  #%>% #split apart version, for better profiling
  pkg_message(c("i"="Making fabletools model"))
    fbm <-fabletools::model(.data=dfmerge, arima = fable::ARIMA(degree),
                      naive = fable::NAIVE(degree))
    fc <- fabletools::forecast(object=fbm, h = h)
    fit <- full_join(fc, dffreq, by="vertex")

  # Fit ARIMA for total edges
  pkg_message(c("i"="Running fabletools on total number of edges"))
  #this fabletools/ARIMA code seems to fail to fit sometimes
  #fit_total <-  tibble::as_tibble(dfall_sum) %>%
  #  tsibble::as_tsibble(index = time) %>%
  #  fabletools::model(arima = fable::ARIMA(edges)) %>% # ~ time
  #  fabletools::forecast(h = h)
  #so use forecast/auto.arima for now
  fit_total <- forecast::auto.arima(dfall_sum$edges) %>% forecast::forecast(h = h, level = conf_level2)

  # Get hilo for vertices separately
  dfhilo <- fit %>%
    fabletools::hilo(level = conf_level2)
  uniq_times <- unique(dfhilo$time)

  colnames(dfhilo)[7] <- 'hilow'
  dfhilo <- dfhilo %>%
    mutate(lower = floor(hilow$lower), upper = ceiling(hilow$upper))
  dfhilo_vtx <- dfhilo %>%
    filter(time == uniq_times[h]) %>%
    group_by(vertex) %>%
    summarize(lower2 = max(min(lower, na.rm = TRUE), 0),
              upper2 = max(upper, na.rm = TRUE),
              mean2 = round(mean(abs(.data$.mean), na.rm = TRUE)),
              freq = mean(freq))

  # Get hilo for total edges
  #pkg_message("Getting hilo of total edges")
  #dfhilo_total <- fit_total %>%
  #  fabletools::hilo(level = conf_level2)


  # Forecast of the old nodes
  # upper limit
  f_rhs_up <- dfhilo_vtx %>%
    pull(upper2)
  # lower limit
  f_rhs_lo <- dfhilo_vtx %>%
    pull(lower2)
  # mean
  f_rhs_mean <- dfhilo_vtx %>%
    pull(mean2)
  # frequency
  freq <- dfhilo_vtx %>%
    pull(freq)
  # vertex label
  vertex <- dfhilo_vtx %>%
    pull(vertex)

  inf_inds <- which(is.infinite(f_rhs_up))
  if(length(inf_inds) > 0){
    f_rhs_up[inf_inds] <- 3*f_rhs_mean[inf_inds]
  }

  #  Total edges predictions
  #total_edges_mean <- dfhilo_total %>%
  #  pull(.data$.mean) %>%
  #  ceiling()
  #colnames(dfhilo_total)[5] <- 'hilow'
  #dfhilo_total <- dfhilo_total %>%
  #  mutate(lower = floor(hilow$lower), upper = ceiling(hilow$upper))
  #total_edges_lower <- dfhilo_total %>%
  #  pull(lower)
  #total_edges_upper <- dfhilo_total %>%
  #  pull(upper)
  #forecast::forecast() code
  total_edges_lower <- as.integer(fit_total$lower) #forecast gives these as numeric/float values, so we need to convert them to integers
  total_edges_upper <- as.integer(fit_total$upper)
  total_edges_mean <- as.integer(fit_total$mean)


  list(
    vertex = vertex,
    degree_hi = f_rhs_up,
    degree_lo = f_rhs_lo,
    degree_mean = f_rhs_mean,
    freq = freq,
    total_edges_mean = total_edges_mean[h],
    total_edges_lower = total_edges_lower[h],
    total_edges_upper = total_edges_upper[h]

  )

}


predict_new_nodes_degree <- function(graphlist){
  # Forecast the degree of new nodes
  # For new nodes forecast the average degree of new nodes

  index <- avg_deg <- forecast <- NULL
  NN <- length(graphlist)
  average_new_degree <-node_count <- edge_count <- rep(0, NN)
  for(jj in 1:NN){
    node_count[jj] <- igraph::vcount(graphlist[[jj]])
    edge_count[jj] <- igraph::ecount(graphlist[[jj]])
    new_node_count <- diff(node_count)
    if(jj == 1){
      average_new_degree[jj] <- 0
    }else{
      average_new_degree[jj] <- mean(igraph::degree(graphlist[[jj]])[(node_count[jj] - new_node_count[(jj-1)] + 1):node_count[jj]])
    }
  }
  df2 <- data.frame(index = 1:NN, avg_deg = average_new_degree)

  new_deg <- df2 %>% tsibble::as_tsibble(index = index) %>%
    fabletools::model(arima = fable::ARIMA(avg_deg)) %>%
    fabletools::forecast(h = 1) %>%
    mutate(forecast = ceiling(.data$.mean)) %>%
    pull(forecast)
  new_deg
}

get_new_nodes_degrees <- function(graphlist, h, confidence_level){
  #get the degree that new nodes should have, based on their age
  num_graphs <- length(graphlist)
  node_counts <- sapply(graphlist, igraph::vcount)
  num_new_nodes <- pmax(diff(node_counts), 0) #number of new nodes each iteration of the graphlist
  new_nodes <- t(mapply(function(graph, nn) c((igraph::vcount(graph)-nn+1),igraph::vcount(graph)), graph=graphlist[-1], nn=num_new_nodes))

  total_entries <- ((num_graphs - 1)*(num_graphs))/2
  new_edge_counts <- matrix(data=0, nrow=total_entries, ncol=3)
  colnames(new_edge_counts) <- c("time", "age", "count")
  out_pos <- 1
  for (i in 2:num_graphs){
    degrees <- igraph::degree(graphlist[[i]])
    average_degrees <- apply(new_nodes[1:i-1,,drop=FALSE], 1,
                              FUN=function(inds){
                                mean(degrees[inds[1]:inds[2]])
                              })
    this_new_entries <- length(average_degrees)
    out_idxs <- out_pos:(out_pos+this_new_entries-1)
    new_edge_counts[out_idxs,1] <- i
    new_edge_counts[out_idxs,2] <- this_new_entries:1
    new_edge_counts[out_idxs,3] <- average_degrees
    out_pos <- out_pos + this_new_entries
  }
  #new nodes edge counts as a tsibble dataframe
  new_node_ec_ts <- as.data.frame(new_edge_counts) %>% tsibble::as_tsibble(key=age, index=time)
  mdl <- fabletools::model(new_node_ec_ts, arima=fable::ARIMA(count))
  fc <- fabletools::forecast(mdl, h=h)
  hilo <- fabletools::hilo(fc, level=confidence_level)
  this_timestep <- filter(hilo, time==(num_graphs + h))
  colnames(this_timestep)[NCOL(this_timestep)] <- "conf"
  degree_map <- matrix(data=0, nrow=NROW(this_timestep), ncol=2)
  degree_map[,1] <- this_timestep$age
  degree_map[,2] <- ceiling(this_timestep$conf$upper) #need to use $ to get the data, [] doesn't work
  degree_map
}

new_nodes_edge_constraints <- function(new_nodes_per_timestep, new_nodes_degrees_per_age){
  #get the degree constraints for new nodes, based on when they should be added to the predicted graph
  #new_nodes_per_timestep is the number of new nodes that should be added at each time step, up to the predicted one
  #new_nodes_degrees_per_age is a 2-column matrix of (age, degree for age) pairs
  total_new_nodes <- sum(new_nodes_per_timestep)
  current_time_step <- length(new_nodes_per_timestep)
  new_nodes_degrees <- rep(0,total_new_nodes)
  out_pos <- 1
  for (i in 1:current_time_step){
    age <- min(current_time_step - i + 1, NROW(new_nodes_degrees_per_age))
    degree_for_age <- new_nodes_degrees_per_age[new_nodes_degrees_per_age[, 1] == (age), 2]
    nodes_added_this_timestep <- new_nodes_per_timestep[[i]]
    new_nodes_degrees[out_pos:(out_pos+nodes_added_this_timestep)] <- degree_for_age
    out_pos <- out_pos + nodes_added_this_timestep
  }
  new_nodes_degrees
}

run_lpsolve <- function(f_con, f_obj, f_rhs, signs){
  f_dir <- signs
  lpobj <- lpSolve::lp("max", f_obj, f_con, f_dir, f_rhs, all.bin = TRUE)
  lpobj
}

run_lpsolve_dense <- function(f_con, f_obj, f_rhs, signs){
  f_con <- as.matrix(f_con)
  f_dir <- signs
  lpobj <- lpSolve::lp(direction = "max",
                       objective.in = f_obj,
                       const.dir = f_dir,
                       const.rhs = f_rhs,
                       dense.const = f_con,
                       all.bin = TRUE)
  lpobj
}

growth_metrics <- function(graphlist){
  n <- length(graphlist)
  deg_cosine_sim <- rep(0, (n-1))
  edges <- vertices <- rep(0, n)
  for( i in 1:n){
    vertices[i] <- igraph::vcount(graphlist[[i]])
    edges[[i]] <- igraph::ecount(graphlist[[i]])
  }
  for(i in 1:(n-1)){
    deg_i <- igraph::degree(graphlist[[i]])
    deg_ip1 <- igraph::degree(graphlist[[i+1]])[1:length(deg_i)]
    deg_cosine_sim[i] <- (sum(deg_i*deg_ip1))/(norm(as.matrix(deg_i), type = "F")*norm(as.matrix(deg_ip1), type = "F"))

  }
  vertex_growth <- vertices[2:n]/vertices[1:(n-1)] - 1
  edges_growth <- edges[2:n]/edges[1:(n-1)] - 1

  list(
    v_growth = vertex_growth,
    e_growth = edges_growth,
    v_growth_mean = mean(vertex_growth[11:(n-1)]),
    v_growth_sd = sd(vertex_growth[11:(n-1)]),
    e_growth_mean = mean(edges_growth[11:(n-1)]),
    e_growth_sd = sd(edges_growth[11:(n-1)]),
    deg_cosine_sim = deg_cosine_sim
  )
}

eval_metrics <- function(actualgr, predgr){

  new_vert <- NULL
  metric1 <- (igraph::vcount(predgr) - igraph::vcount(actualgr))/igraph::vcount(actualgr)

  if(metric1 > 0){
    # vcount(predgr) > vcount(actualgr)
    new_vert <- igraph::vcount(predgr) - igraph::vcount(actualgr)
    actualgr <-igraph::add_vertices(actualgr,nv =new_vert)
  }else if(metric1 < 0){
    # vcount(predgr) < vcount(actualgr)
    new_vert <-  igraph::vcount(actualgr) - igraph::vcount(predgr)
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

sparse_solver_formulation <- function(union_graph, max_degrees, total_edges_constraint){
  #create the LP Solve inputs from a union graph. Each edge in this graph should be considered as a possible component in the predicted graph
  degrees <- igraph::degree(union_graph)
  used_nodes <- which(degrees > 0)
  num_zero_edge_nodes <- igraph::vcount(union_graph) - length(used_nodes)
  pkg_debug(c("i"=sprintf("The union graph has %d nodes, %d nodes have no edges", igraph::vcount(union_graph), num_zero_edge_nodes)))
  objective_function <- igraph::E(union_graph)$weight#igraph::V(union_graph)[used_nodes]
  constraint_rhs <- c(max_degrees[used_nodes], total_edges_constraint)
  num_nodes <- length(used_nodes)
  num_edges <- igraph::ecount(union_graph)
  edge_list <- igraph::as_edgelist(union_graph)
  #remove nodes/constraints with zero edges because LPSolve does not accept them
  if (num_zero_edge_nodes > 0){
    vertex_offset <- 1:(igraph::vcount(union_graph)) - cumsum(degrees <= 0) #offset the from and to vertices by this amount, to remove nodes with zero edges
    edge_list[,1] <-vertex_offset[edge_list[,1]]
    edge_list[,2] <-vertex_offset[edge_list[,2]]
  }

  constraint_matrix_rows <- 3*num_edges
  constraint_matrix <- matrix(data=1, nrow=constraint_matrix_rows, ncol=3)
  #set degree constraints for the source vertex in each edge
  idxs <- 1:num_edges
  constraint_matrix[idxs,1] <- edge_list[,1]
  constraint_matrix[idxs,2] <- idxs
  #set degree constraints for the target vertex in each edge
  end_idxs <- (num_edges+1):(2*num_edges)
  constraint_matrix[end_idxs,1] <- edge_list[,2]
  constraint_matrix[end_idxs,2] <- idxs
  #set total number of edges constraint
  end_idxs <- (2*num_edges+1):(3*num_edges)
  constraint_matrix[end_idxs,1] <- length(constraint_rhs)
  constraint_matrix[end_idxs,2] <- idxs
  constraint_direction = c(rep("<=", num_nodes), "<=") #degree constraints, and total number of edges constraint

  stopifnot(all(constraint_matrix[,1] <= igraph::vcount(union_graph)+1) )
  stopifnot(all(constraint_matrix[,2] <= igraph::ecount(union_graph)  ) )

  list(
    objective_function = objective_function,
    constraint_matrix = constraint_matrix,
    constraint_direction = constraint_direction,
    constraint_rhs = constraint_rhs,
    used_edges = igraph::as_edgelist(union_graph)#use original vertex indices
  )
}

#construct a graph from the LPSolve Solution
graph_from_solution <- function(solution, edge_list, num_nodes){
  indices <- which(solution > 0)
  edges_in_solution <- edge_list[indices,]
  solution_graph <- igraph::graph_from_edgelist(edges_in_solution, directed=FALSE)
  #add in the required number of nodes, in case the last node(s) have zero edges
  if (num_nodes > igraph::vcount(solution_graph)){
    pkg_debug(c("i"=sprintf("adding %d new nodes to solution", num_nodes - igraph::vcount(solution_graph)) ))
    solution_graph <- igraph::add_vertices(solution_graph, num_nodes - igraph::vcount(solution_graph))
  }
  pkg_debug(c("i"=sprintf("The solution graph has %d nodes and %d edges", igraph::vcount(solution_graph), igraph::ecount(solution_graph))))
  solution_graph
}

add_new_nodes_and_edges <- function(graph, new_nodes_degrees, existing_degree_constraints){
  #add the new nodes, and the edges that should be attached to the new nodes
  spare_degrees <- existing_degree_constraints - igraph::degree(graph)
  num_nodes <- igraph::vcount(graph)
  num_new_nodes <- length(new_nodes_degrees)
  graph <- igraph::add_vertices(graph, num_new_nodes)
  for (i in 1:length(new_nodes_degrees)){
    this_node <- num_nodes + i
    this_connections <- sample(num_nodes, new_nodes_degrees[[i]], prob=spare_degrees)
    new_edges <- c(rbind(rep(this_node, new_nodes_degrees[[i]]), this_connections))
    graph <- igraph::add_edges(graph, new_edges)
    spare_degrees[this_connections] <- spare_degrees[this_connections] - 1
  }
  graph
}

check_constraints <-function(graph, degree_constraints, total_edges_constraint){
  #check whether or not a predicted graph meets the given constraints - uses the rhs of the constraints
  num_edges <- length(degree_constraints)
  graph_degrees <- igraph::degree(graph)
  correct <- TRUE
  if (length(graph_degrees) != length(degree_constraints)){
    pkg_message(sprintf("x"="The graph has %d nodes, but there are %d degree constraints", length(graph_degrees), length(degree_constraints)))
    correct <- FALSE
  }

  same_degrees <- graph_degrees <= degree_constraints
  if (!all(same_degrees)) {
    correct <- FALSE
    for (i in 1:num_edges){
      if (graph_degrees[i] > degree_constraints[i]){
        pkg_message(c=("x"=sprintf("Node %d has degree %d but should be <= %d", i, graph_degrees[i], degree_constraints[i])))
      }
    }
  }
  actual_num_edges <- igraph::ecount(graph)
  if (total_edges_constraint != actual_num_edges){
    correct <- FALSE
    pkg_message(c("x"=sprintf("The graph has %d edges, was expecting %d edges", actual_num_edges, total_edges_constraint)))
  }
  if (correct){
    pkg_debug(c("i"="All constraints met") )
  }
  return(correct)
}

#print messages, if in verbose or debug mode
pkg_message <- function(...){
  verbosity <- getOption("netseer.verbose", "quiet")
  if (verbosity %in% c("verbose", "debug")){
    rlang::local_options(rlib_message_verbosity = "verbose")
    rlang::inform(...)
  }
}

#print messages, if in debug mode
pkg_debug <- function(...){
  verbosity <- getOption("netseer.verbose", "quiet")
  if (verbosity == "debug"){
    rlang::local_options(rlib_message_verbosity = "verbose")
    rlang::inform(...)
  }
}
