struct System{N}
    matrix_entries::SparseMatrixCSC{Entry, Int64}
    vector_entries::Vector{Entry}
    diagonal_inverses::Vector{Entry}
    acyclic_children::Vector{Vector{Int64}} # Contains direct children that are not part of a cycle
    cyclic_children::Vector{Vector{Int64}}  # Contains direct and indirect children that are part of a cycle
    parents::Vector{Vector{Int64}}          # Contains direct and cycle-opening parents
    dfs_list::SVector{N,Int64}
    graph::SimpleGraph{Int64}
    dfs_graph::SimpleDiGraph{Int64}
    dimrow::Vector{Int64}
    dimcol::Vector{Int64}

    function System{T}(A, dimrow, dimcol; force_static = false) where T
        N = length(dimrow)
        @assert N == length(dimcol)

        static = force_static || (all(dimrow.<=10) && all(dimcol.<=10))
        full_graph = Graph(A)
        matrix_entries = spzeros(Entry,N,N)

        for i = 1:N
            for j = 1:N
                if i == j
                    matrix_entries[i,j] = Entry{T}(dimrow[i], dimcol[j], static = static)
                elseif j ∈ all_neighbors(full_graph,i)
                    matrix_entries[i,j] = Entry{T}(dimrow[i], dimcol[j], static = static)
                    matrix_entries[j,i] = Entry{T}(dimrow[j], dimcol[i], static = static)
                end
            end
        end

        vector_entries = [Entry{T}(dim, static = static) for dim in dimrow];
        # this is only well-defined for dimrow == dimcol
        diagonal_inverses = [Entry{T}(dim, dim, static = static) for dim in dimrow];

        graphs, roots = split_adjacency(A)
        dfs_list = Int64[]
        acyclic_children = [Int64[] for i=1:N]
        cycles = [Vector{Int64}[] for i=1:N]
        parents = [Int64[] for i=1:N]
        edgelist = LightGraphs.SimpleEdge{Int64}[]

        for (i,graph) in enumerate(graphs)
            root = roots[i]
            dfs_graph = dfs_tree(graph, root)
            sub_dfs_list, cycle_closures = dfs(graph, root)
            append!(dfs_list, sub_dfs_list)

            cycle_dfs_graph = copy(dfs_graph)
            cycle_dfs_graph_reverse = copy(dfs_graph)
            for cycle_closure in cycle_closures
                add_edge!(cycle_dfs_graph, cycle_closure...)
                add_edge!(cycle_dfs_graph_reverse, reverse(cycle_closure)...)
            end
            append!(edgelist, collect(edges(cycle_dfs_graph_reverse)))

            for v in sub_dfs_list
                acyclic_children[v] = neighbors(dfs_graph, v)
                parents[v] = neighbors(reverse(dfs_graph), v)
            end

            for cyclic_members in simplecycles(cycle_dfs_graph)
                v, cyclic_children = cycle_parent_children(cyclic_members, parents)
                cyclic_children = reverse(cyclic_children)
                lump_cycles!(cycles[v], cyclic_children)

                acyclic_children[v] = setdiff(acyclic_children[v], cyclic_children)
                for c in cyclic_children
                    matrix_entries[v,c] = Entry{T}(dimrow[v], dimcol[c], static = static);
                    matrix_entries[c,v] = Entry{T}(dimrow[c], dimcol[v], static = static);

                    v ∉ parents[c] && push!(parents[c],v)
                end
            end
        end

        full_dfs_graph = SimpleDiGraph(edgelist)
        cyclic_children = [unique(vcat(cycles[i]...)) for i=1:N]

        new{N}(matrix_entries, vector_entries, diagonal_inverses, acyclic_children, cyclic_children, parents, dfs_list, full_graph, full_dfs_graph, dimrow, dimcol);
    end

    System(A, dimrow, dimcol; force_static = false) = System{Float64}(A, dimrow, dimcol; force_static = force_static);
end

@inline children(system, v) = outneighbors(system.dfs_graph, v)
@inline connections(system, v) = neighbors(system.graph, v)
@inline parents(system, v) = inneighbors(system.dfs_graph, v)

# There probably exists a smarter way of getting the dense matrix from the spares one
function full_matrix(system::System{N}) where N
    dimrow = system.dimrow
    dimcol = system.dimcol

    range_row = [1:dimrow[1]]
    for (i,dim) in enumerate(collect(Iterators.rest(dimrow, 2)))
        push!(range_row,sum(dimrow[1:i])+1:sum(dimrow[1:i])+dim)
    end
    range_col = [1:dimcol[1]]
    for (i,dim) in enumerate(collect(Iterators.rest(dimcol, 2)))
        push!(range_col,sum(dimcol[1:i])+1:sum(dimcol[1:i])+dim)
    end

    A = zeros(sum(dimrow),sum(dimcol))

    for (i,row) in enumerate(system.matrix_entries.rowval)
        col = findfirst(x->i<x,system.matrix_entries.colptr)-1
            A[range_row[row],range_col[col]] = system.matrix_entries[row,col].value
    end
    return A
end

full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
