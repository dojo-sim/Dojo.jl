# constraints
""" 
    get_joint(mechanism, name) 

    return JointConstraint from Mechanism 

    mechanism: Mechanism 
    name: unique identifier for joint
"""
function get_joint(mechanism::Mechanism, name::Symbol)
    for joint in mechanism.joints
        if joint.name == name
            return joint
        end
    end
    return
end

get_joint(mechanism::Mechanism, id::Integer) = mechanism.joints[id]


# bodies
""" 
    get_body(mechanism, name) 

    returns Body from Mechanism

    mechanism: Mechanism 
    name: unique identifier for body
"""
function get_body(mechanism::Mechanism, name::Symbol)
    if name == :origin
        return mechanism.origin
    else
        for body in mechanism.bodies
            if body.name == name
                return body
            end
        end
    end
    return
end

get_body(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = id == 0 ? mechanism.origin : mechanism.bodies[id-Ne]

# contacts
""" 
    get_contact(mechanism, name) 

    returns ContactConstraint from Mechanism 

    mechanism: Mechanism 
    name: unique identifier for contact
"""
function get_contact(mechanism::Mechanism, name::Symbol)
    for contact in mechanism.contacts
        if contact.name == name
            return contact
        end
    end
    return
end

get_contact(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = mechanism.contacts[id-Ne-Nb]


# nodes
"""
    get_node(mechanism, name) 

    return Node from Mechanism 

    mechanism: Mechanism 
    name: unique identifier for node 
"""
function get_node(mechanism::Mechanism, name::Symbol)
    node = get_body(mechanism, name)
    if node === nothing
        node = get_joint(mechanism,name)
    end
    if node === nothing
        node = get_contact(mechanism,name)
    end
    return node
end

function get_node(mechanism::Mechanism{T,Nn,Ne,Nb}, id::Integer; 
    origin::Bool=false) where {T,Nn,Ne,Nb}
    (origin && id == 0) && return mechanism.origin
    if id <= Ne
        return get_joint(mechanism, id)
    elseif id <= Ne+Nb
        return get_body(mechanism, id)
    else
        return get_contact(mechanism, id)
    end
end

# maximal 
"""
    get_maximal_state(mechanism) 

    return the current maximal state of mechanism 

    mechanism: Mechanism 
"""
function get_maximal_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		x2 = body.state.x2
		v15 = body.state.v15
		q2 = body.state.q2
		ω15 = body.state.ω15
		pack_maximal_state!(z, x2, v15, q2, ω15, i)
	end
	return z
end

"""
    get_next_state(mechanism) 

    return the maximal state of mechanism after one simulation step

    mechanism: Mechanism 
"""
function get_next_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	z_next = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
        x3, v25, q3, ω25 = next_configuration_velocity(body.state, timestep)
		z_next[13 * (i-1) .+ (1:13)] = [x3; v25; vector(q3); ω25]
	end
	return z_next
end

# minimal 
"""
    get_minimal_state(mechanism) 

    return minimal state for mechanism 

    mechanism: Mechanism
"""
function get_minimal_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	x = []

	mechanism = deepcopy(mechanism)
	timestep = mechanism.timestep

    for joint in mechanism.joints
        c = zeros(T,0)
		v = zeros(T,0)
		pbody = get_body(mechanism, joint.parent_id)
		cbody = get_body(mechanism, joint.child_id)
		for element in (joint.translational, joint.rotational)
			pos = minimal_coordinates(element, pbody, cbody)
			vel = minimal_velocities(element, pbody, cbody, timestep)
			push!(c, pos...)
			push!(v, vel...)
		end
		push!(x, [c; v]...)
    end
	x = [x...]
    
	return x
end

function get_minimal_coordinates(mechanism::Mechanism{T}) where T
    d = Dict{Int,Vector{T}}()
    for joint in mechanism.joints
        push!(d, joint.id => minimal_coordinates(mechanism, joint))
    end
    return d
end
