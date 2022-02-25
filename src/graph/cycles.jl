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
    for (i, cycle) in enumerate(cycles)
        if cyclic_children ⊆ cycle
            return
        elseif cyclic_children ⊇ cycle
            append!(inds, i)
        end
    end
    push!(cycles, cyclic_children)
    deleteat!(cycles, inds)
    return
end
