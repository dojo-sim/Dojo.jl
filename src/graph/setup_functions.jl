function dfs(graph,v)
    n = nv(graph)
    visited = zeros(Bool, n)
    list = Int64[]
    cycleclosures = Vector{Int64}[]
    _dfs(graph,v,0,list,cycleclosures,visited)

    # removes double entries of cycle connections and keeps the first found pair
    cycleclosures = cycleclosures[sortperm(sort.(cycleclosures))[1:2:length(cycleclosures)]]

    return list, cycleclosures
end

function _dfs(graph,v,p,list,cycleclosures,visited)
    if !visited[v]
        visited[v] = true
        for node in all_neighbors(graph,v)
            if node != p && visited[node]
                push!(cycleclosures,[v;node])
            end
            _dfs(graph,node,v,list,cycleclosures,visited)
        end
        push!(list,v)
    end

    return
end

function split_adjacency(A)
    graphs = SimpleGraph{Int64}[]

    subinds = connected_components(Graph(A))
    for subset in subinds
        subA = zeros(Int64,size(A)...)
        subA[subset,subset] = A[subset,subset]
        push!(graphs, Graph(subA))
    end

    roots = [subinds[i][1] for i=1:length(subinds)]

    return graphs, roots
end

function cycle_parent_children(cyclic_members, parents)
    parent = -1
    for member in cyclic_members
        !isempty(parents[member]) && parents[member][1] ∈ cyclic_members && continue
        parent = member
        break
    end

    return parent, setdiff(cyclic_members, parent)
end

function lump_cycles!(cycles, cyclic_children)
    inds = Int64[]
    for (i,cycle) in enumerate(cycles)
        if cyclic_children ⊆ cycle
            return
        elseif cyclic_children ⊇ cycle
            append!(inds,i)
        end
    end
    push!(cycles, cyclic_children)
    deleteat!(cycles, inds)
    return
end
