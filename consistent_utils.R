suppressWarnings(library(igraph))
suppressWarnings(library(dequer))

SimpleGraph <- setRefClass(
    "SimpleGraph",
    fields = list(
        order = "numeric",
        adj_matrix_oriented = "matrix",
        adj_matrix = "matrix",
        consistent_candidates = "array",
        #matrix_consistent_candidates = "array",
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
            if(is.null(adj_matrix_)) {
                stop("not implemented")
            }
            else{
                if(is.numeric(adj_matrix_)) adj_matrix_ <- adj_matrix_ > 0
                if(oriented){
                    adj_matrix_oriented <<- adj_matrix_
                    adj_matrix <<- adj_matrix_oriented | t(adj_matrix_oriented)
                }
                else{
                    adj_matrix <<- adj_matrix_
                }
                order <<- ncol(adj_matrix)
            }
            #initialize_consistent_candidates()
            #initialize_matrix_consistent_candidates()
            graph <<- graph_from_adjacency_matrix(adj_matrix, mode="undirected")
            size <<- sum(adj_matrix) / 2
            bcc()
        },

        bcc_aux = function(u, st) {
            # Count of children in current node
            children = 0

            # Initialize discovery time and low value
            Time <<- Time + 1
            depth[u] <<- Time
            lowest[u] <<- Time

            # Recur for all the vertices adjacent to this vertex
            for(v in which(adj_matrix[u,]==TRUE)){
                # If v is not visited yet, then make it a child of u
                # in DFS tree and recur for it
                if(depth[v] == -1){
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
                    if((parent[u] == -1 && children > 1) || (parent[u] != -1 && lowest[v] >= depth[u])){
                        count <<- count + 1 # increment count
                        if(!is_cp[u]){
                            cp_list <<- c(cp_list, u)
                            is_cp[u] <<- 1
                        }

                        w = -1
                        bcc_set = c()
                        while(any(w != c(u,v))){
                            w = pop(st)
                            bcc_set = union(bcc_set, w)
                            #print(w, end='-'),
                        }
                        bcc_list[[length(bcc_list)+1]] <<- bcc_set
                        #print("\n------")
                    }
                }

                else if(v != parent[u] && depth[u] > depth[v]){
                    # Update low value of 'u' only if 'v' is still in stack
                    # (i.e. it's a back edge, not forward edge).
                    # Case 2
                    # -- per Strongly Connected Components Article

                    lowest[u] <<- min(lowest[u], depth[v])
                    push(st, c(u,v))
                }
            }
        },


        bcc = function(){
            Time     <<- 0
            count    <<- 0
            bcc_list <<- list()
            cp_list  <<- c()
            is_cp    <<- rep(0, order)

            depth   <<- rep(-1, order)
            lowest  <<- rep(-1, order)
            parent  <<- rep(-1, order)
            st = stack() # from dequer. modified when passed to function

            for(i in 1:order){
                if(depth[i] == -1)
                    bcc_aux(i, st)

                # If stack is not empty, pop all edges from stack
                if (length(st)>0){
                    count <<- count + 1
                    bcc_set = c()

                    while(length(st)>0){
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
            for(index in 1:length(bcc_list)){
                bcc_set = bcc_list[[index]]
                bc_tree_index = bc_tree_index + 1
                rep = bc_tree_index
                bc_tree_inverse_index[rep] <<- index

                for(node in bcc_set){
                    bcc_set_indices[[node]] <<- c(bcc_set_indices[[node]], index)

                    if(is_cp[node]){
                        if(bc_tree_rep[node] == -1){
                            bc_tree_index = bc_tree_index + 1
                            bc_tree_rep[node] <<- bc_tree_index
                            bc_tree_node_is_cp[bc_tree_index] <<- 1
                            bc_tree_inverse_index[bc_tree_index] <<- node
                        }

                        bc_tree_adj_list[[bc_tree_rep[node]]] <<- c(bc_tree_adj_list[[bc_tree_rep[node]]], rep)
                        bc_tree_adj_list[[rep]] <<- c(bc_tree_adj_list[[rep]], bc_tree_rep[node])

                    }
                    else{
                        bc_tree_rep[node] <<- rep
                    }
                }
            }
        },


        get_candidate_z = function(x, y) {
            if(degree(x) < 1 || degree(y) < 1)
                return(c())

            bcc_common = base::intersect(bcc_set_indices[[x]], bcc_set_indices[[y]])

            if(length(bcc_common)>0){
                return_list = bcc_list[[bcc_common[1]]]
                return_list = return_list[return_list != x & return_list != y]
                #bcc_common = bcc_common[2:length(bcc_common)]
                #return(return_list)
            }
            else{
                start = bc_tree_rep[x]
                end   = bc_tree_rep[y]
                paths = c()
                for(s in bc_tree_bfs(start, end)){
                    if(bc_tree_node_is_cp[s]==0)
                        paths = c(paths, bcc_list[[bc_tree_inverse_index[s]]])
                }
                return_list = setdiff(unique(paths), c(x,y))
            }
            neighbors = base::union(which(adj_matrix[x,]==TRUE), which(adj_matrix[y,]==TRUE))
            return(base::intersect(return_list, neighbors))
        },


        bc_tree_bfs = function(start, end){

            queue = deque()
            pqueue = deque()
            pushback(queue, start)
            pushback(pqueue, c(start))
            visited = c()
            while(length(queue)>=1){
                node = pop(queue)
                path = pop(pqueue)
                visited = union(visited, node)
                for(neighbor in base::setdiff(bc_tree_adj_list[[node]], c(path, visited)) ) {
                    if(neighbor == end) {
                        return(c(path, neighbor))
                    }
                    else{
                        pushback(queue, neighbor)
                        pushback(pqueue, c(path, neighbor))
                    }
                }
            }
        },


        bfs = function(start, end, excludes=NULL) {
            if(start %in% excludes || end %in% excludes) {
                return(NULL)
            }
            queue = deque()
            pqueue = deque()
            pushback(queue, start)
            pushback(pqueue, c(start))
            visited = c()
            while(length(queue)>=1){
                node = pop(queue)
                path = pop(pqueue)
                visited = union(visited, node)
                for(neighbor in base::setdiff(neighbors(graph, node),
                                            c(path, excludes, visited)) ) {
                    if(neighbor == end) {
                        return(c(path, neighbor))
                    }
                    else{
                        pushback(queue, neighbor)
                        pushback(pqueue, c(path, neighbor))
                    }
                }
            }
        },

        bfs_collider_free = function(start, end, excludes=NULL, has_one_way=FALSE) {
            if(start %in% excludes || end %in% excludes) {
                return(NULL)
            }
            queue = deque()
            pqueue = deque()
            pushback(queue, start)
            pushback(pqueue, c(start))
            visited = c()
            while(length(queue)>=1){
                has_one_way = FALSE
                node = pop(queue)
                path = pop(pqueue)
                if(length(path) > 5) return(NULL)
                if((length(path) < 2) || ! is_one_way(path[length(path)-1], path[length(path)])){
                    visited = c(visited, node)
                    if(length(path)>=2){
                    has_one_way = FALSE
                    }
                }
                else
                    has_one_way = TRUE
                candidate_neighbors = which(adj_matrix_oriented[node,])
                if(! has_one_way)
                    candidate_neighbors = base::union(candidate_neighbors, which(adj_matrix_oriented[,node]==1))
                candidate_neighbors = base::setdiff(
                    candidate_neighbors,
                    base::union(path, base::union(excludes, visited)))
                for(neighbor in candidate_neighbors) {
                    if(neighbor == end) {
                        return(c(path, neighbor))
                    }
                    else{
                        pushback(queue, neighbor)
                        pushback(pqueue, c(path, neighbor))
                    }
                }
            }
        },

        is_one_way = function(head, tail){
            return(adj_matrix_oriented[head,tail]==1 && adj_matrix_oriented[tail,head]==0)
        },

        dfs_potential_collider_free = function(start, end, path=NULL, visited=c(),
                                     excludes=c(), has_one_way=FALSE){
            #"""path from start to end without collider"""

            if(is.null(path))
                path = c(start)
            if(start == end)
                return(path)
            if((length(path) < 2) || ! is_one_way(path[length(path)-1], path[length(path)]))
                visited = c(visited, start)
            else
                has_one_way = TRUE
            candidate_neighbors = which(adj_matrix_oriented[start,])
            if(! has_one_way)
                candidate_neighbors = base::union(candidate_neighbors, which(adj_matrix_oriented[,start]==1))
            candidate_neighbors = base::setdiff(
                candidate_neighbors,
                base::union(path, base::union(excludes, visited)))

            for(neighbor in candidate_neighbors){
                res = dfs_potential_collider_free(neighbor, end, c(path,neighbor),
                                                  visited, excludes, has_one_way)
                if(!is.null(res)) return(res)
            }
        },

        dfs_collider_free = function(start, end, path=NULL, visited=c(),
                                     excludes=c(), has_one_way=FALSE){
            #"""path from start to end without collider"""

            if(is.null(path))
                path = c(start)
            if(start == end)
                return(path)
            if( (length(path) < 2) || ! is_one_way(path[length(path)-1], path[length(path)]))
                visited = c(visited, start)
            else
                has_one_way = TRUE
            candidate_neighbors = which(adj_matrix_oriented[start,])
            if(! has_one_way)
                candidate_neighbors = base::union(candidate_neighbors, which(adj_matrix_oriented[,start]==1))
            candidate_neighbors = base::setdiff(
                candidate_neighbors,
                base::union(path, base::union(excludes, visited)))

            for(neighbor in candidate_neighbors){
                res = dfs_collider_free(neighbor, end, c(path,neighbor),
                                        visited, excludes)
                if(!is.null(res)) return(res)
            }
        },

        dfs_one_way = function(start, end, path=NULL, visited=c(),
                               excludes=c()){

            if(is.null(path))
                path = c(start)
            if(start == end)
                return(path)
            visited = c(visited, start)
            candidate_neighbors = which(adj_matrix_oriented[start,])
            #candidate_neighbors = base::union(candidate_neighbors, which(adj_matrix[,start]==1))
            candidate_neighbors = base::setdiff(
                candidate_neighbors,
                base::union(path, base::union(excludes, visited)))

            for(neighbor in candidate_neighbors){
                res = dfs_one_way(neighbor, end, c(path,neighbor),
                                  visited, excludes)
                if(!is.null(res)) return(res)
            }
        },

        degree = function(x){
            return(sum(adj_matrix[x,]))
        },

        initialize_consistent_candidates = function(){
            consistent_candidates <<- array(-1, ncol(adj_matrix))
        },

        #initialize_matrix_consistent_candidates = function(){
        #    matrix_consistent_candidates <<- array(rep(-1, ncol(adj_matrix)^3), dim=rep(ncol(adj_matrix),3))
        #},

        is_consistent = function(x, y, z, oriented=TRUE, debug=FALSE, check_potential_colliders=FALSE){

            if(length(z)>1){
                for(single_z in z){
                    if(!is_consistent(x,y,single_z,debug=debug)) return(FALSE)
                }
                return(TRUE)
                #return(all(sapply(z, function(z){is_consistent(x, y, z)})))
            }
            if (x == z || y == z
                || degree(z) <= 1
                || degree(x) <= 0
                || degree(y) <= 0) {
                #consistent_candidates[z] <<- 0
                #matrix_consistent_candidates[x,y,z] <<- 0
                if(debug) print("1")
                return(FALSE)
            }

            #if(consistent_candidates[z] != -1) {
            #    return(consistent_candidates[z]==1)
            #}
            #if(matrix_consistent_candidates[x,y,z] != -1) {
            #    return(matrix_consistent_candidates[x,y,z]==1)
            #}

            is_consistent = FALSE

            if(!oriented){
                return(z %in% get_candidate_z(x,y))
            }
            if(check_potential_colliders){
                if(adj_matrix[x,z]){
                    if(is_one_way(x,z)){
                        return(!is.null(dfs_one_way(z,y,excludes=c(x))))
                    }
                    else{
                        return(!is.null(dfs_potential_collider_free(z,y,excludes=c(x))))
                    }
                }
                else{
                    if(!adj_matrix[y,z]){
                        # Z is not a neighbor of X or Y
                        return(FALSE)
                    }
                    if(is_one_way(y,z)){
                        return(!is.null(dfs_one_way(z,x,excludes=c(y))))
                    }
                    else{
                        return(!is.null(dfs_potential_collider_free(z,x,excludes=c(y))))
                    }
                }
            }
            else{
                if(!(z %in% get_candidate_z(x,y))){
                    if(debug) print("2")
                    return(FALSE)
                }
                else{
                    # if z is not a neighbor of neither x or y, return false
                    if(!adj_matrix[x,z] && !adj_matrix[y,z]) {
                        if(debug) print("3")
                        return(FALSE)
                    }
                    # if z is a non-child of x or y
                    if (adj_matrix_oriented[z,y]==1 || adj_matrix_oriented[z,x]==1){
                        return(TRUE)
                    }
                    else{
                        if(debug) print("4")
                        return(FALSE)
                    }
                    #if(adj_matrix[x,z]){
                        #if(adj_matrix[y,z]){
                        #    if(!(is_one_way(x,z) && is_one_way(y,z))) return(TRUE)
                        #}

                        #is_consistent <- (!is.null(bfs_collider_free(z,y,excludes=c(x),has_one_way=is_one_way(x,z))))
                    #}
                    #else{
                    #    if(!adj_matrix[y,z]){
                    #        # Z is not a neighbor of X or Y
                    #        is_consistent <- (FALSE)
                    #    }
                        #is_consistent <- (!is.null(bfs_collider_free(z,x,excludes=c(y),has_one_way=is_one_way(y,z))))
                    #}
                }
            }
            #consistent_candidates[z] <<- is_consistent
            #matrix_consistent_candidates[x,y,z] <<- is_consistent
            #return(is_consistent)
        },


        add_non_oriented_edge = function(x, y){
            adj_matrix[x,y] <<- TRUE
            adj_matrix[y,x] <<- TRUE
            adj_matrix_oriented[x,y] <<- TRUE
            adj_matrix_oriented[y,x] <<- TRUE
        }


        #dfs_collider_free_2 = function(start, end, excludes=NULL) {
        #    if(start %in% excludes || end %in% excludes) {
        #        return(NULL)
        #    }
        #    queue = deque()
        #    pqueue = deque()
        #    pushback(queue, start)
        #    pushback(pqueue, c(start))
        #    visited = c()
        #    while(length(queue)>=1){
        #        node = popback(queue)
        #        path = popback(pqueue)
        #        visited = union(visited, node)
        #        for(neighbor in base::setdiff(neighbors(graph, node),
        #                                    c(path, excludes, visited)) ) {
        #            if(neighbor == end) {
        #                return(c(path, neighbor))
        #            }
        #            else{
        #                pushback(queue, neighbor)
        #                pushback(pqueue, c(path, neighbor))
        #            }
        #        }
        #    }
        #}

    )
)
