function velocity_index(mechanism::Mechanism{T,Nn,Ne}) where {T,Nn,Ne}
    ind = []
    off = 0
	for id in mechanism.root_to_leaves
        (id > Ne) && continue # only treat joints
        joint = mechanism.joints[id]
        nu = control_dimension(joint)
        push!(ind, Vector(off + nu .+ (1:nu)))
        off += 2nu
    end
    return vcat(ind...)
end

# find all the joints parents of a body
function parent_joints(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, body::Body) where {T,Nn,Ne,Nb,Ni}
	ids = parents(mechanism.system, body.id)
	ids = intersect(ids, 1:Ne) # filter out the bodies
	return [get_node(mechanism, id) for id in ids]
end

