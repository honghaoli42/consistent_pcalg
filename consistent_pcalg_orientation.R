setClass("pcAlgoConsistent", contains = "pcAlgo",
         slots = c(GAfterInitialization = "matrix")) ## zMin for compatibility


skeleton <- function(suffStat, indepTest, alpha, labels, p,
		     method = c("stable", "original", "stable.consistent"), m.max = Inf,
		     fixedGaps = NULL, fixedEdges = NULL,
		     NAdelete = TRUE, numCores = 1, verbose = FALSE,
             lastState = NULL, initial_ord=0L)
{
  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of DAG given data
  ## Order-independent version! NEU
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat: List containing all necessary elements for the conditional
  ##             independence decisions in the function "indepTest".
  ## - indepTest: predefined function for testing conditional independence
  ## - alpha: Significance level of individual partial correlation tests
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - numCores: number of cores to be used for calculation if
  ##   method = "stable.fast"
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G, sepset, pMax, ord, n.edgetests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009
  ## Modification: Diego Colombo; Martin Maechler; Alain Hauser

  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p))
      p <- length(labels)
    else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    ## Don't want message, in case this is called e.g. from fciPlus():
    ## else
    ##   message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
                                        ##  if (method == "stable.fast" && .Platform$OS.type == "windows") {
                                        ##    method <- "stable"
                                        ##    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
                                        ##  }

  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps

  diag(G) <- FALSE

  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")

  GAfterInitialization = matrix()
  ## Check number of cores
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && numCores > 0)
  if (numCores > 1 && method != "stable.consistent") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  if (method == "stable.consistent") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose),
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
        NAdelete = NAdelete,
        numCores = numCores)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    ## sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
  }
  else {
    ## Original R version

    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- initial_ord
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,"\n",sep = "")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remEdges) {
        if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        #if(!is.null(lastState)) lastState$initialize_consistent_candidates()
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
          nbrsBool[y] <- FALSE
          if(is.null(lastState)) nbrs <- seq_p[nbrsBool]
          else{
            consistent_Zs = lastState$get_candidate_z(x,y)
            nbrs <- base::intersect(seq_p[nbrsBool], consistent_Zs) #Get only consistent candidates for separation set
            non_descendants <- which(lastState$adj_matrix_oriented[,x]==1)
            nbrs <- base::intersect(nbrs, non_descendants)
          }
          #nbrs <- sample(seq_p[nbrsBool], length(nbrs), replace=F) #Get random candidate set the size of consistent set
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              S_is_consistent = TRUE
              #if((!is.null(lastState)) && ord>0){
              #  if(pval >= alpha) S_is_consistent = lastState$is_consistent(x,y,S)
              #  #if(!S_is_consistent)print(paste(x,y, paste(S, collapse=",")))
              #}
              if((pval >= alpha && S_is_consistent)) { # independent and consistent with orientation with relation to last state
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      if(ord==0) GAfterInitialization = G
      ord <- ord + 1L
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }

  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }

  ## final object
  new("pcAlgoConsistent", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1),
      GAfterInitialization = GAfterInitialization)
}## end{ skeleton }



pc <- function(suffStat, indepTest, alpha, labels, p,
               fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
               u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.consistent"),
               conservative = FALSE, maj.rule = FALSE,
               solve.confl = FALSE, numCores = 1, verbose = FALSE)
{
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - gTrue: Graph suffStatect of true DAG
  ## - conservative: If TRUE, conservative PC is done
  ## - numCores: handed to skeleton(), used for parallelization
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Diego Colombo, Markus Kalisch

  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")

    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }

  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")


  ################################################################################
  ### CONSISTENCY

  no_loop_graph_states = TRUE
  graphObject = NULL
  first_iteration = TRUE
  GAfterInitialization = NULL
  graph_states = list()
  consistent_iter = 0

  while(no_loop_graph_states && consistent_iter < 100){

    ## Skeleton
    skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                    fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                    NAdelete=NAdelete, m.max=m.max, numCores=numCores, verbose=verbose,
                    lastState=graphObject, initial_ord = ifelse(first_iteration, 0, 1))

    # Save adjacency matrix with unconditional (in)dependences to use as starting point in future iterations
    if(first_iteration) {
        GAfterInitialization = skel@GAfterInitialization
        if(is.null(fixedGaps)){ fixedGaps = !GAfterInitialization}
        else{ fixedGaps = fixedGaps | (!GAfterInitialization)}
    }

    ## Orient edges
    orientation <- NULL
    if (!conservative && !maj.rule) {
        orientation <- switch (u2pd,
                               "rand" = udag2pdag(skel),
                               "retry" = udag2pdagSpecial(skel)$pcObj,
                               "relaxed" = udag2pdagRelaxed(skel, verbose = verbose, solve.confl = solve.confl))
    }
    else { ## u2pd "relaxed" : conservative _or_ maj.rule
        ## version.unf defined per default
        ## Tetrad CPC works with version.unf=c(2,1)
        ## see comment on pc.cons.intern for description of version.unf
        pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                              version.unf = c(2,1), maj.rule = maj.rule, verbose = verbose)
        orientation <- udag2pdagRelaxed(pc.$sk, verbose = verbose,
                                        unfVect = pc.$unfTripl, solve.confl = solve.confl)
    }

    G = as(orientation@graph, "matrix")
    graphObject = SimpleGraph(G, oriented=TRUE)

    sepSets = orientation@sepset

    # Save current consistent graph state
    G = graphObject$adj_matrix_oriented
    print("# EDGES :")
    print( sum((G+t(G))>0)/2)

    # Check if the obtained adjacency matrix G is a previous result
    if(!first_iteration){
        for(statei in 1:length(graph_states)) {
            if(all(graph_states[[statei]] == G)){
                # Limit cycle found
                print("Loop detected")
                print(statei)
                no_loop_graph_states = F
                for(statek in statei:length(graph_states)){
                    # Get the union of all adjacency matrices in the limit cycle
                    G = G | graph_states[[statek]]
                    #sepset becomes wrong
                }
                break
            }
        }
    }

    graph_states[[length(graph_states)+1]] = G
    first_iteration = FALSE
    #is_graph_consistent = graphObject$is_graph_consistent()
    consistent_iter = consistent_iter + 1
  }
  print("# EDGES FINAL:")
  print( sum((G+t(G))>0)/2)
  graphObject = SimpleGraph(G, oriented=TRUE)

  #skel@call <- cl # so that makes it into result
  ################################################################################

  # Resolve sepset after consistency loop union
  for(x in 1:(p-1)){
    for(y in (x+1):p){
      if(G[x,y] || G[y,x]){
        sepSets[[x]][[y]] <- numeric(0)
        sepSets[[y]][[x]] <- numeric(0)
      }
    }
  }

  sepSets_matrix = do.call(rbind, sepSets)
  which_sepSets = which(matrix(vapply(sepSets_matrix, function(x){!(is.null(x) ||
                                                             length(x)==0)}, logical(1)),
                               ncol=p), arr.ind=T)

  #print(orientation@sepset)
  #print(vapply(sepSets_matrix, function(x){!is.null(x) || length(x)!=0}, logical(1)))
  #print(which_sepSets)
  #print(sepSets[[1]][[30]])
  # Check consistency of all separation sets in the union of the limit cycle
  if(nrow(which_sepSets)>0){
    add_back_list = c()
    for(sepSet_index in 1:nrow(which_sepSets)){

        #if(as(skel@graph, 'matrix')[x,y]) stop("Edge exists with sepset")
        sepSet = which_sepSets[sepSet_index,]
        x = sepSet[1]
        y = sepSet[2]
        z = sepSets[[x]][[y]]
        # Checks if Z is consistent for X,Y pair, adds back X-Y edge if not
        if(!graphObject$is_consistent(x,y,z)){
            add_back_list = c(add_back_list, sepSet_index) # List of inconsistent pairs with respect to current graph
            print(G[x,y] || G[y,x])
            print(z)
            print(paste(x,y,paste(z,collapse=',')))
            graphObject$is_consistent(x,y,z,debug=TRUE)
        }
    }

    # Add back inconsistent pairs all at once
    print(paste("Adding back", length(add_back_list), "edges."))
    if(length(add_back_list)>0){
        write("Non zero add_back_list !", "~/Documents/PhD/Work/consistent_simulation/alert.txt", append=TRUE)
    }
    for(sepSet_index in 1:length(add_back_list)){
        sepSet = which_sepSets[add_back_list[sepSet_index],]
        x = sepSet[1]
        y = sepSet[2]
        graphObject$add_non_oriented_edge(x,y)
    }
    graphObject$bcc()
    #graphObject$initialize_matrix_consistent_candidates()
  }


  result = orientation
  result@graph <- as(graphObject$adj_matrix_oriented, "graphNEL")
  result@sepset <- sepSets

  result

} ## {pc}