function global_id!(nodes) 
    oldnewid = Dict([node.id => i for (i, node) in enumerate(nodes)]...)
    for node in nodes
        node.id = oldnewid[node.id]
        if typeof(node) <: Constraint
            node.parent_id = get(oldnewid, node.parent_id, 0)
            node.child_id = get(oldnewid, node.child_id, 0)
        end
    end
end