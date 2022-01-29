function global_id!(nodes) 
    oldnewid = Dict([node.id => i for (i, node) in enumerate(nodes)]...)
    for node in nodes
        node.id = oldnewid[node.id]
        if typeof(node) <: Constraint
            node.parentid = get(oldnewid, node.parentid, nothing)
            node.childids = [get(oldnewid, childid, nothing) for childid in node.childids]
        end
    end
end