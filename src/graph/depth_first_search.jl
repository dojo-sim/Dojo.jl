function dfs(graph, v)
    n = nv(graph)
    visited = zeros(Bool, n)
    list = Int64[]
    cycleclosures = Vector{Int64}[]
    dfs!(graph, v, 0, list, cycleclosures, visited)

    # removes double entries of cycle connections and keeps the first found pair
    cycleclosures = cycleclosures[sortperm(sort.(cycleclosures))[1:2:length(cycleclosures)]]

    return list, cycleclosures
end

function dfs!(graph, v, p, list, cycleclosures, visited)
    if !visited[v]
        visited[v] = true
        for node in all_neighbors(graph, v)
            if node != p && visited[node]
                push!(cycleclosures, [v; node])
            end
            dfs!(graph, node, v, list, cycleclosures, visited)
        end
        push!(list, v)
    end

    return
end