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
                                   dense_opt = 2,
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

  debug <- FALSE#TRUE#
  print(sprintf("Using weights option %d", weights_opt))
  #print(sprintf("Starting RHS (max degrees):"))
  #print(probj$degree_hi)
  print(sprintf("Max degrees: %d", probj$total_edges_upper))
  lpsolve_vars_type <- 2
  if (debug | lpsolve_vars_type == 1){#TRUE
    # Step 1 - forecast the number of nodes - done in predict_graph function


    # Step 2 - Forecast the degree of the old nodes - done in predict_graph function
    # But the other hilo values need to be taken from probj
    f_rhs_hi <- probj$degree_hi
    f_rhs_lo <- probj$degree_lo
    freq <- probj$freq
    vertex <- probj$vertex
    total_hi <- probj$total_edges_upper
    total_low <- probj$total_edges_lower

    # Step 3 - Forecast the degree of new nodes
    # For new nodes forecast the average degree of new nodes
    rm_verts <- new_deg <- NULL
    rhsobj <- compute_rhs(graphlist, formulation, new_nodes, f_rhs_hi, f_rhs_lo = NULL, vertex, freq, total_hi, total_low)
    f_rhs <- rhsobj$f_rhs
    new_deg <- rhsobj$new_deg
    rm_verts <- rhsobj$rm_verts
    signs <- rhsobj$signs
    check_rhs <- f_rhs


    # Step 4 - Construct the union graph, find the constraint matrix and the weights
    biggr <- construct_union_graph(graphlist, new_nodes, new_deg, rm_verts, weights_opt = weights_opt, weights_param = weights_param)
    if(dense_opt == 1){
      f_con <- constraint_matrix(biggr)
    }else if(dense_opt == 2){
      f_con <- constraint_matrix_dense(biggr)
    }
    print(sprintf("f_con first has %d unique nodes, max %d, has dim %d", length(unique(f_con[,1])), max(f_con[,1]), dim(f_con)))
    wts <- get_weights(biggr, f_con, weights_opt, dense_opt)
    #print(sprintf("wts len %d", length(wts)) )
    #print(wts)
    # Step 4.1 - Remove zero weights and associated entries for optimization
    tighter_obj <- remove_zero_rows_and_columns(wts, f_rhs, f_con, f_obj, signs, dense_opt)
    f_con <- tighter_obj$f_con
    f_rhs <- tighter_obj$f_rhs
    signs <- tighter_obj$signs
    f_obj <- tighter_obj$f_obj
    print(sprintf("f_con later has %d unique nodes, max %d, dim %s", length(unique(f_con[,1])), max(f_con[,1]), dim(f_con)))
    print(sprintf("f_obj len %d", length(f_obj)))
    #print(f_obj)


    if(dense_opt == 1){
      # Step 5 - Remove zero rows and Run lpSolve
      # Remove zero rows
      rhs_nonzero <- remove_zero_rhs(f_con, f_rhs, signs)
      f_con <- rhs_nonzero$f_con
      f_rhs <- rhs_nonzero$f_rhs
      signs <- rhs_nonzero$signs
      # # Remove >= rows. Because it didn't work for one example. NOT HAPPY ABOUT THIS!
      # rhs_geq <- remove_less_than(f_con, f_rhs, signs)
      # f_con <- rhs_geq$f_con
      # f_rhs <- rhs_geq$f_rhs
      # signs <- rhs_geq$signs
    }
    print(sprintf("orig method has %d nodes, %d potential edges, %d edge constraint", igraph::vcount(biggr), length(f_obj), f_rhs[length(f_rhs)]))
    #print(sprintf("orig weights"))
    #print(f_obj)
    obj_orig <- f_obj
    #print("orig weights before new nodes")
    #print(f_obj[length(f_obj)-new_nodes:length(f_obj)])
    #print("orig rhs:")
    #print(f_rhs)
    #print("orig signs/dir:")
    #print(signs)
    #print("orig constraitns:")
    #print(f_con)

    if(dense_opt == 1){
      lpobj <-  run_lpsolve(f_con, f_obj, f_rhs, signs)
    }else if(dense_opt == 2){
      lpobj <-  run_lpsolve_dense(f_con, f_obj, f_rhs, signs)
    }


    # Step 6: Get the graph relating to the values
    # Using the solution get the bit string - add the zeros that were there originally
    bitstring <- rep(0, length(wts))
    bitstring[which(wts!=0)] <- lpobj$solution
    #bitstring <- lpobj$solution #testing all edges

    if(dense_opt == 1){
      adj <- adjacency_from_solution(bitstring, igraph::gorder(biggr))
      grout <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    }else if (dense_opt == 2){
      edges <- edgelist_from_solution(bitstring, igraph::gorder(biggr))
      grout <- igraph::graph_from_edgelist(edges, directed = FALSE)
      add_verts <- igraph::vcount(biggr) - igraph::vcount(grout)
      if(add_verts > 0){
        grout <- igraph::add_vertices(grout, add_verts)
      }


    }
    print(sprintf("orig method produced a graph with %d nodes and %d edges", igraph::vcount(grout), igraph::ecount(grout)))

    plot(grout)
    title("OriginAl Graph Setup")
    gr_orig <- grout
    lpobj_orig <- lpobj
  }
  if (debug | lpsolve_vars_type == 2){ #TRUE else#
    print("Setting maximum degree constraint")
    max_degrees <- probj$degree_hi
    if (new_nodes > 0){
      new_degrees <- predict_new_nodes_degree(graphlist)
      maximum_degrees <- c(max_degrees, rep(new_degrees, new_nodes)) #the maximum degrees for each vertex in the predicted graph
    } else {
      maximum_degrees <- max_degrees
    }
    #print("Max degrees:")
    #print(maximum_degrees)
    print("Edges in each graph:")
    print(unlist(lapply(graphlist, igraph::ecount)))
    print("Constructing union graph")
    biggr <- construct_union_graph(graphlist, new_nodes, maximum_degrees, rm_verts, weights_opt = weights_opt, weights_param = weights_param)
    print(sprintf("The union graph has %d nodes, %d edges, constraint to use %d edges", igraph::vcount(biggr), igraph::ecount(biggr), probj$total_edges_upper))
    print("Setting sparse solver constraints")
    lpsolve_vars <- sparse_solver_formulation(biggr, maximum_degrees, probj$total_edges_upper)#probj$degree_hi
    #print(lpsolve_vars)
    #print(lpsolve_vars$constraint_rhs)
    #print(lpsolve_vars$objective_function)
    obj_new <- lpsolve_vars$objective_function
    l <- length(lpsolve_vars$objective_function)
    print("Running LPSolve")
    lpobj <- run_lpsolve_dense(lpsolve_vars$constraint_matrix, lpsolve_vars$objective_function, lpsolve_vars$constraint_rhs, lpsolve_vars$constraint_direction)
    print(sprintf("LPSolve solution status (should be 0) is %d", lpobj$status) )
    #print(lpobj$solution)
    #print(lpobj)
    print("Constructing graph from LPSove Solution")
    grout <- graph_from_solution(lpobj$solution, lpsolve_vars$used_edges, igraph::vcount(biggr))
    check_rhs <- lpsolve_vars$constraint_rhs
    #plot(grout)
    #title("Alternete Graph Setup")
    gr_new <- grout
  }
  print("Checking that constraints are met")
  check_constraints(grout, maximum_degrees, probj$total_edges_upper)

  if (debug){
    #print(gr_orig)
    #print(gr_new)
    old_new <- igraph::difference(gr_orig, gr_new)
    #plot(old_new)
    #title(sprintf("old - new: %d edges", igraph::ecount(old_new)) )
    print(sprintf("old - new: %d edges", igraph::ecount(old_new)))
    new_old <- igraph::difference(gr_new, gr_orig)
    #plot(new_old)
    #title(sprintf("new - old: %d edges", igraph::ecount(new_old)) )
    print(sprintf("new - old: %d edges", igraph::ecount(new_old)) )
    print(eval_metrics(gr_orig, gr_new))
    print(sprintf("orig obj len %d, new obj function length %d", length(obj_orig), length(obj_new)))
    #print("objective function new - objective function original:")
    #print(obj_new - obj_orig)
    print(sprintf("Orig objective: %f, new objective %f: difference = %f", lpobj_orig$objval, lpobj$objval, lpobj$objval - lpobj_orig$objval))
  }

  print("Predict graph internal done")
  grout

}

remove_zero_rows_and_columns <- function(wts, f_rhs, f_con, f_obj, signs, dense_opt){

  V1 <- V2 <- V3 <- renamed <- original <- NULL
  #f_obj <- wts
  # Columns - keep only non-zero ones
  nonzero_col_inds <- which(wts != 0)
  edge_index_lookup <- 1:length(wts) - cumsum(wts <= 0)
  f_obj <- wts[nonzero_col_inds]
  use_constraints <- f_con[,2] %in% nonzero_col_inds
  f_con <- f_con[use_constraints,]
  f_con[,2] <- edge_index_lookup[f_con[,2]]
  col_renamed_df <- data.frame(original = nonzero_col_inds, renamed = 1:length(nonzero_col_inds))

  # Rows - keep only non-zero ones
  nonzero_row_inds <- which(f_rhs !=0)
  print(sprintf("rhs has len %d, nonzero len %d", length(f_rhs), length(nonzero_row_inds)))
  vertex_index_lookup <- 1:(length(f_rhs) ) - cumsum(f_rhs <= 0)
  f_rhs <- f_rhs[nonzero_row_inds]
  signs <- signs[nonzero_row_inds]
  f_con[,1] <- vertex_index_lookup[f_con[,1]]
  row_renamed_df <- data.frame(original = nonzero_row_inds, renamed = 1:length(nonzero_row_inds))
  print("row renaming:")
  print(row_renamed_df)

  # Remove and rename the zero columns
  if(dense_opt == 1){
    f_con <- f_con[ ,nonzero_col_inds]
  }else if(dense_opt == 2){
    #f_con <- f_con %>%
    #  filter(V3 !=0)
    #f_con <- f_con %>%
    #  filter(V2 %in% col_renamed_df$original) %>%
    #  rename(original = V2) %>%
    #  left_join(col_renamed_df) %>%
    #  rename(V2 = renamed ) %>%
    #  select(-original) %>%
    #  select(V1, V2, V3)
  }


  # Remove and rename the zero rows
  #if(dense_opt == 1){
  #  f_con <- f_con[nonzero_row_inds, ]
  #}else if(dense_opt == 2){
  #  f_con <- f_con %>%
  #    filter(V1 %in% row_renamed_df$original) %>%
  #    rename(original = V1) %>%
  #    left_join(row_renamed_df) %>%
  #    rename(V1 = renamed) %>%
  #    select(-original) %>%
  #    select(V1, V2, V3)
  #}

  list(
    f_con = f_con,
    f_rhs = f_rhs,
    signs = signs,
    f_obj = f_obj)
}


get_weights <- function(biggr, f_con, weights_opt, dense_opt){
  # The length of the objective function
  weight <- From <- To <- value <- constraint <- column <- NULL

  if(dense_opt == 1){
    len <- dim(f_con)[2]
  }else if(dense_opt == 2){
    len <- max(f_con[ ,2])
  }

  # weights schemes
  f_obj <- rep(0, len)
  if(weights_opt == 1){
    f_obj <- rep(1, len)
  }else if((weights_opt == 2)|(weights_opt == 3)){
    if(dense_opt == 1){
      # standard constraint matrix
      colsums <- apply(f_con, 2, sum)
      f_obj[which(colsums != 0)] <- 1
    }else if(dense_opt == 2){
      # dense constraint matrix
      colnames(f_con) <- c("constraint", "column", "value")
      cols <- f_con %>%
        filter(value > 0) %>%
        distinct(column) %>%
        pull(column)
      f_obj[cols] <- 1
    }
  }else if(weights_opt == 4 | weights_opt == 5 | weights_opt == 6 | weights_opt == 7){
    # Weights df
    weights <- cbind.data.frame(igraph::as_edgelist(biggr), igraph::E(biggr)$weight)
    colnames(weights) <- c("From", "To", "weight")
    weights <- weights %>%
      arrange(From, To)
    # Edges df as in the constraint matrix
    n <- igraph::vcount(biggr)
    from_seq <- c()
    to_seq <- c()
    for( i in 1:(n-1)){
      t_seq <- (i+1):n
      fr_seq <- rep(i, length(t_seq))
      from_seq <- c(from_seq, fr_seq)
      to_seq <- c(to_seq, t_seq)
    }
    df_edges <- data.frame("From" = from_seq, "To" = to_seq)
    f_obj <- left_join(df_edges, weights) %>%
      arrange(From, To) %>%
      tidyr::replace_na(list(weight = 0)) %>%
      pull(weight)
  }
  f_obj
}



remove_zero_rhs <- function(f_con, f_rhs, signs){
  inds <- which(f_rhs == 0)
  if(length(inds) > 0){
    f_rhs <- f_rhs[-inds]
    f_con <- f_con[-inds, ]
    signs <- signs[-inds]
  }
  list(
    f_con = f_con,
    f_rhs = f_rhs,
    signs = signs)

}

remove_less_than <- function(f_con, f_rhs, signs){
  inds <- which(signs == '>=')
  if(length(inds) > 0){
    f_rhs <- f_rhs[-inds]
    f_con <- f_con[-inds, ]
    signs <- signs[-inds]
  }
  list(
    f_con = f_con,
    f_rhs = f_rhs,
    signs = signs)

}

compute_rhs <- function(graphlist,
                        formulation,
                        new_nodes,
                        f_rhs_hi,
                        f_rhs_lo,
                        vertex,
                        freq,
                        total_hi,
                        total_lo){

  degree <- new_deg <- f_rhs <- rm_verts <- NULL

  if(new_nodes > 0){
    new_deg <- predict_new_nodes_degree(graphlist)
    if(is.null(f_rhs_lo)){
      f_rhs <- c(f_rhs_hi, rep(new_deg, new_nodes))
      signs <- c(rep("<=", length(f_rhs_hi)),  rep("<=", new_nodes) )
    }else{
      f_rhs <- c(f_rhs_hi, rep(new_deg, new_nodes), f_rhs_lo, rep(0, new_nodes))
      signs <- c(rep("<=", length(f_rhs_hi)),  rep("<=", new_nodes), rep(">=", length(f_rhs_lo)), rep(">=", new_nodes) )
    }
  }else if(new_nodes == 0){
    f_rhs <- c(f_rhs_hi, f_rhs_lo)
    signs <- c(rep("<=", length(f_rhs_hi)), rep(">=", length(f_rhs_lo)) )
  }else if(new_nodes < 0){
    rm_num <- -1*new_nodes
    signs <- c(rep("<=", length(f_rhs_hi)), rep(">=", length(f_rhs_lo)) )
    df_temp <- data.frame(vertex = vertex, degree = c(f_rhs_hi, f_rhs_lo), signs = signs, freq = freq)
    df_temp_rm <- df_temp %>%
      filter(degree == 0) %>%
      arrange(freq) %>%
      dplyr::slice(1:rm_num)
    df_temp2 <- dplyr::setdiff(df_temp, df_temp_rm)
    f_rhs <- df_temp2 %>%
      pull(degree)
    rm_verts <- df_temp_rm %>%
      pull(vertex)
    signs <- df_temp2 %>%
      pull(signs)
  }

  # If only formulation = 2
  # Add the upper bound for total edges
  if(formulation == 2){
    f_rhs <- c(f_rhs, total_hi)
    signs <- c(signs, "<=")
  }

  list(
    new_deg = new_deg,
    f_rhs = f_rhs,
    rm_verts = rm_verts,
    signs = signs
  )
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
  list(
    new_nodes = new_nodes[h],
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
      print(sprintf("Processing graph %d", ii))
      grr <- graphlist[[ii]]
      igraph::E(grr)$weight <- 1/(NN - ii + 1) # this line make harmonically decaying weights
      if(ii == 1){
        biggr <- grr
      }else{
        biggr <- (biggr %u% grr)
        igraph::E(biggr)$weight_1[is.na(igraph::E(biggr)$weight_1)] <- 0
        igraph::E(biggr)$weight_2[is.na(igraph::E(biggr)$weight_2)] <- 0
        igraph::E(biggr)$weight <- igraph::E(biggr)$weight_1 +  igraph::E(biggr)$weight_2
        #print(igraph::E(biggr)[[1]])
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


  # Add edges if nodes share a neighbour
  print("Adding neighbours of neighbours")
  bigadj1 <- igraph::as_adjacency_matrix(biggr)
  bigadj2 <- bigadj1 %*% bigadj1
  non_edges <- Matrix::which((bigadj1 == 0) & (bigadj2 > 0), arr.ind=TRUE)
  print(sum(non_edges[,1] >= non_edges[,2]))
  non_edges <- t(non_edges[non_edges[,1] < non_edges[,2],])
  #non_l <- list()
  #lcount <- 1
  #for(ll in 1:(NROW(bigadj2)-1)){
  #  for(mm in (ll+1):NCOL(bigadj2)){
  #    if((bigadj2[ll, mm] > 0) & (bigadj1[ll, mm] == 0)){
  #      non_l[[lcount]] <- c(ll, mm)
  #      lcount <- lcount + 1
  #    }
  #  }
  #}
  #non_orig <- matrix(unlist(non_l), byrow=TRUE, ncol=2)
  #nn_order <- order(non_edges[,1], non_edges[,2])
  #non_edges <- non_edges[nn_order,]
  #print("NoN nodes - adj matrix version")
  #print(non_edges)
  #print("Non edges - looping through array version")
  #print(non_orig)
  #print(all(non_edges  ==  non_orig))
  #non_edges <- non_orig
  #non_edges <- unlist(non_l)
  #print("mat and list edges")
  #print(non_edges)
  #print(unlist(non_l))
  #g1 <- igraph::add_edges(biggr, non_edges, weight = 2*new_weights)
  #g2 <- igraph::add_edges(biggr, unlist(non_l), weight = 2*new_weights)
  #plot(igraph::difference(g1, g2))
  #plot(igraph::difference(g2, g1))
  if(weights_opt %in% c(4,5,6,7)){
    biggr <- igraph::add_edges(biggr, non_edges, weight = 2*new_weights)
  }else{
    # weights_opt 1, 2, 3
    biggr <- igraph::add_edges(biggr, non_edges)
  }



  if(new_nodes > 0){
    # Add new nodes
    print("Adding new nodes, and their edges")
    biggr <- igraph::add_vertices(biggr, new_nodes)
    new_vertices <- igraph::V(biggr)[( igraph::gorder(biggr)-new_nodes+1):igraph::gorder(biggr)]
    #old_vertices <- igraph::V(biggr)[1:(igraph::gorder(biggr)-new_nodes)]
    # Use only a fixed number of old vertices
    # Which vertices have the highest degree
    grlast <- graphlist[[NN]]
    num_attach <- 10 # attach potential edges to this number of nodes in the union graph, with the highest degrees
    print(sprintf("grlast has %d nodes and %d edges", igraph::vcount(grlast), igraph::ecount(grlast)))
    if(igraph::vcount(grlast) > num_attach){
      old_vertices <- igraph::V(biggr)[order(igraph::degree(grlast), decreasing = TRUE)[1:num_attach]]
    }else{
      old_vertices <- igraph::V(biggr)[1:(igraph::gorder(biggr)-new_nodes)]
    }

    new_deg2 <- min(2*new_deg, igraph::vcount(biggr))
    verts <- order(igraph::degree(biggr), decreasing = TRUE)[1:new_deg2]

    if((weights_opt == 1)|(weights_opt == 2)){
      # Binary weights - new nodes connected to all old nodes
      # Add new edges from new nodes to all old nodes
      possible_edges <- c(rbind(rep(old_vertices, new_nodes), rep(new_vertices, each = length(old_vertices)) ))
    }else if(weights_opt == 3){
      # Binary weights - new nodes connected to most connected old nodes
      # Add new edges from new nodes to mostly connected vertices
      possible_edges <- c(rbind(rep(verts, new_nodes), rep(new_vertices, each = length(verts)) ))
    }else if(weights_opt == 4|weights_opt == 5|weights_opt == 6|weights_opt == 7){
      # Proportional weights - new nodes connected to all old nodes
      # But the weights will be much smaller
      possible_edges <- c(rbind(rep(old_vertices, new_nodes), rep(new_vertices, each = length(old_vertices)) ))
      # new_weights <- quantile(igraph::E(biggr)$weight, probs = weights_param)
    }

    if(weights_opt == 4|weights_opt == 5|weights_opt == 6|weights_opt == 7){
      biggr <- biggr %>%
        igraph::add_edges(possible_edges,weight = new_weights)
    }else{
      biggr <- biggr %>%
        igraph::add_edges(possible_edges)
    }

  }

  total_edges <- igraph::ecount(biggr)
  new_edges <- total_edges - existing_edges
  print(sprintf("The union graph has %d existing edges, %d new edges, for %d total edges and %.2f%% are new nodes", existing_edges, new_edges, total_edges, 100*new_edges/total_edges))

  biggr
}


predict_old_nodes_degree <- function(graphlist, conf_level2, h){
  # Forecast the degree of the old nodes
  future::plan(future::multisession)
  degree <- edges <- hilow <- lower <- upper <- upper2 <- lower2  <- mean2 <- NULL

  NN <- length(graphlist)
  print("Setting up dfall")
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

  print("Making dfmerge")
  dfmerge <- tibble::as_tibble(dfall) %>%
    tsibble::as_tsibble(index = time , key = vertex) %>%
    dplyr::arrange(time, vertex) %>%
    tsibble::fill_gaps(degree = 0)

  print("making dffreq")
  dffreq <- dfall %>%
    group_by(vertex) %>%
    summarize(freq = n())


  # Fit ARIMA models
  #fit <- dfmerge %>%
  #  fabletools::model(arima = fable::ARIMA(degree),
  #                    naive = fable::NAIVE(degree)) %>%
  #  fabletools::forecast(h = h) %>%
  #  full_join(dffreq)

  #%>% #split apart version, for better profiling
  print("making fabletools model")
    fbm <-fabletools::model(.data=dfmerge, arima = fable::ARIMA(degree), # ~ time
                      naive = fable::NAIVE(degree)) # ~ time
    print("doing fabletools forecast")
    fc <- fabletools::forecast(object=fbm, h = h)#%>%
    print("full join of fabletools results")
    fit <- full_join(fc, dffreq)

  #original version - may handle profiling better
  #fit <- dfmerge %>%
  #  fabletools::model(arima = fable::ARIMA(degree),
  #                    naive = fable::NAIVE(degree)) %>%
  #  fabletools::forecast(h = h) %>%
  #  full_join(dffreq)

  # Fit ARIMA for total edges
    print("running fabletools on total number of edges")
    #print("input:")
    #print(dfall_sum)
    #print("tibble:")
    #print(tibble::as_tibble(dfall_sum))
    #print("model:")
    #print(tibble::as_tibble(dfall_sum) %>%
    #        tsibble::as_tsibble(index = time) %>%
    #        fabletools::model(arima = fable::ARIMA(edges ~ time)))
  #fit_total <-  tibble::as_tibble(dfall_sum) %>%
  #  tsibble::as_tsibble(index = time) %>%
  #  fabletools::model(arima = fable::ARIMA(edges)) %>% # ~ time
  #  fabletools::forecast(h = h)
  fit_total <- forecast::auto.arima(dfall_sum$edges) %>% forecast::forecast(h = h, level = conf_level2)
  print("total edges fit:")
  print(fit_total)

  # Get hilo for vertices separately
  print("getting hilo of vertices")
  dfhilo <- fit %>%
    fabletools::hilo(level = conf_level2)
  uniq_times <- unique(dfhilo$time)

  print("hilo mutate")
  colnames(dfhilo)[7] <- 'hilow'
  dfhilo <- dfhilo %>%
    mutate(lower = floor(hilow$lower), upper = ceiling(hilow$upper))
  print("hilo filter/groupby/summarise")
  dfhilo_vtx <- dfhilo %>%
    filter(time == uniq_times[h]) %>%
    group_by(vertex) %>%
    summarize(lower2 = max(min(lower, na.rm = TRUE), 0),
              upper2 = max(upper, na.rm = TRUE),
              mean2 = round(mean(abs(.data$.mean), na.rm = TRUE)),
              freq = mean(freq))

  # Get hilo for total edges
  #print("Getting hilo of total edges")
  #dfhilo_total <- fit_total %>%
  #  fabletools::hilo(level = conf_level2)


  # Forecast of the old nodes
  print("Getting hilo forecast of old nodes")
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
  print("Getting total edges hilo")
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
  #forecast::forecast code
  total_edges_lower <- as.integer(fit_total$lower) #forecast gives these as numeric/float values
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


constraint_matrix <- function(gr){
  # gr is an igraph
  adj <- igraph::as_adjacency_matrix(gr)
  num_cols <- dim(adj)[1]*(dim(adj)[1] - 1)/2
  NN <- num_rows <- dim(adj)[1]
  constr_mat <- matrix(0, nrow = num_rows, ncol = num_cols)
  st_col <- 1
  for(jj in 1:(num_rows-1)){
    vals <- adj[jj,(jj+1):NN]
    en_col <- st_col + length(vals) - 1
    # Update relevant row in the constraint matrix
    constr_mat[jj, st_col:en_col] <- vals
    # Fill the diagonal in the submatrix below
    if(dim(diag(vals))[1] == 0){
      # Last entry - no matrix, only 1 entry
      constr_mat[(jj+1):NN,st_col:en_col] <- vals
    }else{
      constr_mat[(jj+1):NN,st_col:en_col] <-  diag(vals)
    }

    # Update st_col
    st_col <- en_col + 1

  }

  # All edges constraint
  all_edges <- rep(1, dim(constr_mat)[2])
  constr_mat <- rbind.data.frame(constr_mat, all_edges)

  constr_mat
}


constraint_matrix_dense <- function(gr){
  # gr is an igraph
  From <- To <- V1 <- V2 <- NULL

  edges <- igraph::as_edgelist(gr)
  colnames(edges) <- c("From", "To")
  n <- igraph::vcount(gr)
  n_seq <- c(0,seq(n-1, 1))
  n_seq2 <- cumsum(n_seq)
  edges_1d <- edgelist_to_indices(edges, n)
  #number of rows = dense_mat is sum(max(2*edges from each vertex, 1)) + number of vertices (maybe max of destination vertcex count, there might be vertices that do not have any edges going to them)
  out_rows <- 2*NROW(edges) + (n - length(unique(edges[,1]))) + max(edges_1d)
  dense_matrix <-matrix(data=1, nrow = out_rows, ncol=3) #initialise the matrix to 1
  this_start <- 0
  for( i in 1:n){
    vertex <- i
    #dat_edges <- edges %>% #this is being done for each n, and is fairly expensive
    #  as.data.frame() %>%
    #  dplyr::filter(From == i)
    #this should be faster
    use_edges <- edges[,1] == i #this seems to be much faster than filter
    dat_edges <- as.data.frame(edges[use_edges,,drop=FALSE]) #I don't think this conversion is necessary, but will need to check
    if(NROW(dat_edges) > 0){
      #df <- matrix(0, nrow = 2*NROW(dat_edges), ncol = 3)
      to_edges <- dat_edges[ ,2]
      cols <- n_seq2[i] + to_edges-i #this seems to be the 1d index of the edges
      len <- NROW(dat_edges)
      idxs <- this_start + 1:len
      # The first part of the constraints
      # en <- st + len - 1
      #df[1:len, 1] <- i
      #df[1:len, 2] <- cols
      #df[1:len, 3] <- 1
      #df[(len+1):(2*len), 1] <- to_edges
      #df[(len+1):(2*len), 2] <- cols
      #df[(len+1):(2*len), 3] <- 1
      dense_matrix[idxs, 1] <- i
      dense_matrix[idxs, 2] <- cols
      #dense_matrix[idxs, 3] <- 1
      idxs <- idxs + len
      dense_matrix[idxs, 1] <- to_edges
      dense_matrix[idxs, 2] <- cols
      #dense_matrix[idxs, 3] <- 1
      this_start <- this_start + 2*len
    }else{
      #df <- matrix(c(i, n_seq2[i], 0), nrow = 1)
      #this seems to be setting a constraint to zero, which shouldn't matter - you could skip this entry
      dense_matrix[this_start,] <- c(i, n_seq2[i], 0)
      this_start <- this_start + 1
    }
    #if(i == 1){
    #  dense_mat <- df
    #}else{
    #  dense_mat <- rbind.data.frame(dense_mat, df) #rbinding togther each loop iteration
    #}
  }
  # Total edges constraint
  #doing max(dense_matrix[,2])can return a value from the case for a node with no edges, not sure if this is intended
  num_rows <- max(edges_1d)#max(dense_matrix[ ,2]) #this might be different, this is the cols variable #this should probably be the maximum edge id, not the maximum node id
  #edges_const <- data.frame(constraint = (n+1), coef = 1:num_rows, val = 1)
  coeffs <- 1:num_rows
  idxs <- this_start + coeffs
  dense_matrix[idxs,1] <- n+1
  dense_matrix[idxs,2] <- coeffs
  #dense_matrix[idxs, 3] <- 1
  #colnames(edges_const) <- colnames(dense_mat)
  #dense_mat <- rbind.data.frame(dense_mat, edges_const)

  dense_mat <- dense_matrix %>%
    as.data.frame() %>% #need to check if this conversion to a dataframe is necessary
    arrange(V1, V2) #v1 and v2 do not appear to be set - these appear to be column names which can be used even if they are not set as variables

  dense_mat

}

sparse_solver_formulation <- function(union_graph, max_degrees, total_edges_constraint){
  #create the LP Solve inputs from a union graph. Assuming that each edge in this graph should be considered as a possible component in the predicted graph
  degrees <- igraph::degree(union_graph)
  used_nodes <- which(degrees > 0)
  num_zero_edge_nodes <- igraph::vcount(union_graph) - length(used_nodes)
  print(sprintf("The union graph has %d nodes, %d nodes have no edges", igraph::vcount(union_graph), num_zero_edge_nodes))
  objective_function <- igraph::E(union_graph)$weight#igraph::V(union_graph)[used_nodes]
  constraint_rhs <- c(max_degrees[used_nodes], total_edges_constraint)
  num_nodes <- length(used_nodes)
  num_edges <- igraph::ecount(union_graph)
  edge_list <- igraph::as_edgelist(union_graph)#
  if (num_zero_edge_nodes > 0){
    vertex_offset <- 1:(igraph::vcount(union_graph)) - cumsum(degrees <= 0) #offset the from and to vertices by this amount, to remove nodes with zero edges
    edge_list[,1] <-vertex_offset[edge_list[,1]]
    edge_list[,2] <-vertex_offset[edge_list[,2]]
  }
  do_sort <- FALSE  #optionally sort edge list to match sevvandi's order
  if (do_sort){
    edge_order <- order(edge_list[,1], edge_list[,2])
    edge_list <- edge_list[edge_order,]
    objective_function <- objective_function[edge_order]
  }

  constraint_matrix_rows <- 3*num_edges
  constraint_matrix <- matrix(data=1, nrow=constraint_matrix_rows, ncol=3)
  #print(sprintf("before, max const idx is %d", max(constraint_matrix[,1])))
  #set degree constraints for the source vertex in each edge
  idxs <- 1:num_edges
  constraint_matrix[idxs,1] <- edge_list[,1]
  constraint_matrix[idxs,2] <- idxs
  #print(sprintf("after to edges, max const idx is %d", max(constraint_matrix[,1])))
  #set degree constraints for the target vertex in each edge
  end_idxs <- (num_edges+1):(2*num_edges)
  constraint_matrix[end_idxs,1] <- edge_list[,2]
  constraint_matrix[end_idxs,2] <- idxs
  #print(sprintf("after from edges, max const idx is %d", max(constraint_matrix[,1])))
  #set total number of edges constraint
  end_idxs <- (2*num_edges+1):(3*num_edges)
  constraint_matrix[end_idxs,1] <- length(constraint_rhs)
  constraint_matrix[end_idxs,2] <- idxs
  constraint_direction = c(rep("<=", num_nodes), "<=")

  #print(sprintf("The max from node is %d and the max to node is %d", max(edge_list[,1]), max(edge_list[,2]) ))
  #print(sprintf("Max of constraint col 1 is %d, the last constraint should be %d", max(constraint_matrix[,1]), length(constraint_rhs) ) )
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


constraint_matrix_with_newnodes <- function(gr, newnodes){
  # gr is an igraph
  # newnodes is the number of new nodes
  adj <- igraph::as_adjacency_matrix(gr)
  num_cols <- dim(adj)[1]*(dim(adj)[1] - 1)/2
  NN <- num_rows <- dim(adj)[1]
  constr_mat <- matrix(0, nrow = num_rows, ncol = num_cols)
  new_positions <- rep(0, num_cols)
  st_col <- 1
  for(jj in 1:(num_rows-1)){
    vals <- adj[jj,(jj+1):NN]
    en_col <- st_col + length(vals) - 1
    # Update relevant row in the constraint matrix
    constr_mat[jj, st_col:en_col] <- vals
    # Fill the diagonal in the submatrix below
    if(dim(diag(vals))[1] == 0){
      # Last entry - no matrix, only 1 entry
      constr_mat[(jj+1):NN,st_col:en_col] <- vals
    }else{
      constr_mat[(jj+1):NN,st_col:en_col] <-  diag(vals)
    }

    # New positions vector
    new_positions[(en_col-newnodes + 1):en_col] <- 1

    # Update st_col
    st_col <- en_col + 1

  }
  list(
    constr_mat = constr_mat,
    new_positions = new_positions)

}



adjacency_from_solution <- function(x, n){
  # x is the solution vector
  # n is the number of vertices
  adj <- matrix(0, nrow = n, ncol = n)
  adj[lower.tri(adj)] <- x
  adj2 <- t(adj)
  adj2[lower.tri(adj2)] <- x
  adj2
}

edgelist_from_solution <- function(x, n){
  # x is the solution vector
  # n is the number of vertices
  #print(x)
  #print(n)
  #v_seq <- 1:n
  version <- 2
  if (version == 1){
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
    edges <-as.matrix(dfall[which(x == 1), ])
    return(edges)
  } else {
    idxs <- which(x == 1)
    offset <- 0 #where the current source vertex starts (it starts one after this one)
    out_idx <- 1 #where to save this index
    n <- (sqrt(8*length(x) +1) +1)/2 #this may have accuracy issues for large lengths
    src <- 1 #source vertex
    dst <- 1 #destination vertex
    this_src <- n - src #the number of 1d indices that start at this source vertex
    out <-matrix(data=0, nrow=length(idxs), ncol=2) #output array
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
    return(out)
  }
    #print(all(edges == out) )
    #return(out)

}

edgelist_to_indices <- function(edgelist, num_vertices){
  #convert an edge list (2-column matrix of (from, to) vertex index pairs) into 1d array indexes into an upper triangle array
  #assumes that the second element is > than the first element
  indices <- (1-edgelist[,1])*(edgelist[,1] - 2*num_vertices)/2 + (edgelist[,2] - edgelist[,1])
  indices
}

graph_from_solution <- function(solution, edge_list, num_nodes){
  #todo - i think i need to implement an integer square root for this to work
  #actually r's sqrt seems to work well for integers
  indices <- which(solution > 0)
  edges_in_solution <- edge_list[indices,]
  solution_graph <- igraph::graph_from_edgelist(edges_in_solution, directed=FALSE)
  if (num_nodes > igraph::vcount(solution_graph)){
    print(sprintf("adding %d nodes", num_nodes - igraph::vcount(solution_graph)) )
    solution_graph <- igraph::add_vertices(solution_graph, num_nodes - igraph::vcount(solution_graph))
  }
  print(sprintf("The solution graph has %d nodes and %d edges", igraph::vcount(solution_graph), igraph::ecount(solution_graph)))
  solution_graph
}

check_constraints <-function(graph, degree_constraints, total_edges_constraint){
  #check whether or not a predicted graph meets the given constraints - uses the rhs of the constraints
  #print(length(constraints))
  #print(constraints)
  #print(length(igraph::degree(graph)))
  #print(igraph::degree(graph))
  #print(igraph::ecount(graph))
  num_edges <- length(degree_constraints)
  # <- constraints[1:num_edges]
  # <- constraints[num_edges+1]
  graph_degrees <- igraph::degree(graph)
  #print(graph_degrees)
  #print(edge_constraints)
  print(sprintf("%d graph degrees, %d degree_constraints", length(graph_degrees), length(degree_constraints)))
  same_degrees <- graph_degrees <= degree_constraints
  correct <- TRUE
  if (!all(same_degrees)) {
    correct <- FALSE
    for (i in 1:num_edges){
      if (graph_degrees[i] > degree_constraints[i]){
        print(sprintf("Node %d has degree %d but should be <= %d", i, graph_degrees[i], degree_constraints[i]))
      }
    }
  }
  actual_num_edges <- igraph::ecount(graph)
  if (total_edges_constraint != actual_num_edges){
    correct <- FALSE
    print(sprintf("The graph has %d edges, was expecting %d edges", actual_num_edges, total_edges_constraint))
  }
  if (correct){
    print("All constraints met")
  }
  return(correct)
}
