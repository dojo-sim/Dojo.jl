function split_adjacency(A)
    graphs = SimpleGraph{Int64}[]

    subinds = connected_components(Graph(A))
    for subset in subinds
        subA = zeros(Int64, size(A)...)
        subA[subset, subset] = A[subset, subset]
        push!(graphs, Graph(subA))
    end

    roots = [subinds[i][1] for i=1:length(subinds)]

    return graphs, roots
end

