function velocity_index(mechanism::Mechanism{T,Nn,Ne}) where {T,Nn,Ne}
    ind = []
    off = 0
	for id in mechanism.root_to_leaves
        (id > Ne) && continue # only treat joints
        joint = mechanism.joints[id]
        nu = input_dimension(joint)
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

# dimensions
function minimal_dimension(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nx = 0
    free_rot_base = false # we are going to check if the link attached to the base has free orientation
    nx = 2 * input_dimension(mechanism, ignore_floating_base = false)
    free_rot_base && (nx += 1)
    return nx
end

maximal_dimension(mechanism::Mechanism{T,Nn,Ne,Nb}; attjac::Bool=false) where {T,Nn,Ne,Nb} = attjac ? 12Nb : 13Nb

function input_dimension(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}; 
    ignore_floating_base::Bool=false) where {T,Nn,Ne,Nb,Ni}
    nu = 0
    for joint in mechanism.joints
        nu += input_dimension(joint, ignore_floating_base=ignore_floating_base)
    end
    return nu
end

