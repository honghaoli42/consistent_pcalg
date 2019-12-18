suppressWarnings(library(igraph))
suppressWarnings(library(dequer))

# Class for biconnected component analysis
SimpleGraph <- setRefClass(
  "SimpleGraph",
  fields = list(
    order = "numeric",
    adj_matrix_oriented = "matrix",
    adj_matrix = "matrix",
    consistent_candidates = "array",
    graph = "ANY",
    size = "numeric",
    Time = "numeric",
    count = "numeric",
    bcc_list = "ANY",
    cp_list = "ANY",
    is_cp = "numeric",
    depth = "numeric",
    lowest = "numeric",
    parent = "numeric",
    bc_tree_adj_list      = "ANY",
    bc_tree_inverse_index = "ANY",
    bc_tree_node_is_cp    = "ANY",
    bcc_set_indices       = "ANY",
    bc_tree_rep           = "ANY"
  ),
  methods = list(

    initialize = function(adj_matrix_, order, oriented=FALSE){
      if (is.null(adj_matrix_)) {
        stop("not implemented")
      }
      else {
        if (is.numeric(adj_matrix_)) adj_matrix_ <- adj_matrix_ > 0
        if (oriented) {
          adj_matrix_oriented <<- adj_matrix_
          adj_matrix <<- adj_matrix_oriented | t(adj_matrix_oriented)
        }
        else {
          adj_matrix <<- adj_matrix_
        }
        order <<- ncol(adj_matrix)
      }
      graph <<- graph_from_adjacency_matrix(adj_matrix, mode="undirected")
      size <<- sum(adj_matrix) / 2
      bcc()
    },

    # Auxiliary recurrent method for biconnected component decomposition.
    bcc_aux = function(u, st) {
      # Count of children in current node
      children = 0
      # Initialize discovery time and low value
      Time <<- Time + 1
      depth[u] <<- Time
      lowest[u] <<- Time
      # Recur for all the vertices adjacent to this vertex
      for (v in which(adj_matrix[u,]==TRUE)) {
        # If v is not visited yet, then make it a child of u
        # in DFS tree and recur for it
        if (depth[v] == -1) {
          parent[v] <<- u
          children = children + 1
          push(st, c(u,v)) # store the edge in stack
          bcc_aux(v, st)
          # Check if the subtree rooted with v has a connection to
          # one of the ancestors of u
          # Case 1 -- per Strongly Connected Components Article
          lowest[u] <<- min(lowest[u], lowest[v])
          # If u is an articulation point, pop
          # all edges from stack till (u, v)
          if ((parent[u] == -1 && children > 1) ||
              (parent[u] != -1 && lowest[v] >= depth[u])) {
            count <<- count + 1 # increment count
            if (!is_cp[u]) {
              cp_list <<- c(cp_list, u)
              is_cp[u] <<- 1
            }
            w = -1
            bcc_set = c()
            while (any(w != c(u,v))) {
              w = pop(st)
              bcc_set = union(bcc_set, w)
            }
            bcc_list[[length(bcc_list)+1]] <<- bcc_set
          }
        }
        else if (v != parent[u] && depth[u] > depth[v]) {
          # Update low value of 'u' only if 'v' is still in stack
          # (i.e. it's a back edge, not forward edge).
          # Case 2
          # -- per Strongly Connected Components Article
          lowest[u] <<- min(lowest[u], depth[v])
          push(st, c(u,v))
        }
      }
    },

    # Biconnected components decomposition of the graph,
    # allowing for search of candidate vertices for separation.
    bcc = function() {
      Time     <<- 0
      count    <<- 0
      bcc_list <<- list()
      cp_list  <<- c()
      is_cp    <<- rep(0, order)
      depth    <<- rep(-1, order)
      lowest   <<- rep(-1, order)
      parent   <<- rep(-1, order)
      st = stack() # from dequer. modified when passed to function

      for (i in 1:order) {
        if (depth[i] == -1)
          bcc_aux(i, st)
        # If stack is not empty, pop all edges from stack
        if (length(st) > 0) {
          count <<- count + 1
          bcc_set = c()
          while (length(st) > 0) {
            w = pop(st)
            bcc_set = union(bcc_set, w)
          }
          bcc_list[[length(bcc_list)+1]] <<- bcc_set
        }
      }

      bc_tree_size = length(cp_list)+length(bcc_list)
      bc_tree_adj_list      <<- vector(mode="list", length=bc_tree_size)
      bc_tree_inverse_index <<- rep(-1, bc_tree_size)
      bc_tree_node_is_cp    <<- rep(0, bc_tree_size)
      bcc_set_indices       <<- vector(mode="list", length=order)
      bc_tree_rep           <<- rep(-1, order)

      bc_tree_index = 0
      for (index in 1:length(bcc_list)) {
        bcc_set = bcc_list[[index]]
        bc_tree_index = bc_tree_index + 1
        rep = bc_tree_index
        bc_tree_inverse_index[rep] <<- index

        for (node in bcc_set) {
          bcc_set_indices[[node]] <<- c(bcc_set_indices[[node]], index)
          if (is_cp[node]) {
            if (bc_tree_rep[node] == -1) {
              bc_tree_index = bc_tree_index + 1
              bc_tree_rep[node] <<- bc_tree_index
              bc_tree_node_is_cp[bc_tree_index] <<- 1
              bc_tree_inverse_index[bc_tree_index] <<- node
            }

            bc_tree_adj_list[[bc_tree_rep[node]]] <<- c(bc_tree_adj_list[[bc_tree_rep[node]]], rep)
            bc_tree_adj_list[[rep]] <<- c(bc_tree_adj_list[[rep]], bc_tree_rep[node])
          }
          else {
            bc_tree_rep[node] <<- rep
          }
        }
      }
    },

    get_candidate_z = function(x, y) {
      if (degree(x) < 1 || degree(y) < 1)
        return(c())
      bcc_common = base::intersect(bcc_set_indices[[x]], bcc_set_indices[[y]])
      if (length(bcc_common) > 0) {
        return_list = bcc_list[[bcc_common[1]]]
        return_list = return_list[return_list != x & return_list != y]
      }
      else {
        start = bc_tree_rep[x]
        end   = bc_tree_rep[y]
        paths = c()
        for (s in bc_tree_bfs(start, end)) {
          if (bc_tree_node_is_cp[s]==0)
            paths = c(paths, bcc_list[[bc_tree_inverse_index[s]]])
        }
        return_list = setdiff(unique(paths), c(x,y))
      }
      neighbors = base::union(which(adj_matrix[x,]==TRUE), which(adj_matrix[y,]==TRUE))
      return(base::intersect(return_list, neighbors))
    },

    # Return the shortest path between two nodes in the block-cut tree.
    bc_tree_bfs = function(start, end) {
      queue = deque()
      pqueue = deque()
      pushback(queue, start)
      pushback(pqueue, c(start))
      visited = c()
      while (length(queue) >= 1) {
        node = pop(queue)
        path = pop(pqueue)
        visited = union(visited, node)
        for (neighbor in base::setdiff(bc_tree_adj_list[[node]], c(path, visited)) ) {
          if(neighbor == end) {
            return(c(path, neighbor))
          }
          else {
            pushback(queue, neighbor)
            pushback(pqueue, c(path, neighbor))
          }
        }
      }
    },

    degree = function(x) {
      return(sum(adj_matrix[x,]))
    },

    is_consistent = function(x, y, z, oriented=TRUE){
      if (length(z) > 1) {
        for (single_z in z) {
          if (!is_consistent(x, y, single_z)) return(FALSE)
        }
        return(TRUE)
      }
      if (x == z || y == z || degree(z) < 2 || degree(x) < 1 || degree(y) < 1) {
        return(FALSE)
      }
      if (!oriented) {
        return(z %in% get_candidate_z(x, y))
      }
      else {
        if (!(z %in% get_candidate_z(x, y))) {
          return(FALSE)
        }
        else {
          # if z is not a neighbor of neither x or y, return false
          if (!adj_matrix[x,z] && !adj_matrix[y,z]) {
            return(FALSE)
          }
          # if z is a non-child of x or y
          if (adj_matrix_oriented[z,y] == 1 || adj_matrix_oriented[z,x] == 1) {
            return(TRUE)
          }
          else {
            return(FALSE)
          }
        }
      }
    },

    add_non_oriented_edge = function(x, y) {
      adj_matrix[x,y] <<- TRUE
      adj_matrix[y,x] <<- TRUE
      adj_matrix_oriented[x,y] <<- TRUE
      adj_matrix_oriented[y,x] <<- TRUE
    }
  )
)
