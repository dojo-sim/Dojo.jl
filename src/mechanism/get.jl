# states
function get_current_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
        x2, v15, q2, ϕ15 = initial_configuration_velocity(body.state)
		z[13 * (i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
	end
	return z
end

function get_next_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep = mechanism.timestep
	z_next = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
        x3, v25, q3, ϕ25 = next_configuration_velocity(body.state, timestep)
		z_next[13 * (i-1) .+ (1:13)] = [x3; v25; vector(q3); ϕ25]
	end
	return z_next
end

# constraints
get_joint_constraint(mechanism::Mechanism, id::Integer) = mechanism.joints[id]

function get_joint_constraint(mechanism::Mechanism, name::Symbol)
    for joint in mechanism.joints
        if joint.name == name
            return joint
        end
    end
    return
end

# bodies
get_body(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = id == 0 ? mechanism.origin : mechanism.bodies[id-Ne]

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

# contacts
get_contact_constraint(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, id::Integer) where {T,Nn,Ne,Nb,Ni} = mechanism.contacts[id-Ne-Nb]

function get_contact_constraint(mechanism::Mechanism, name::Symbol)
    for contact in mechanism.contacts
        if contact.name == name
            return contact
        end
    end
    return
end

# nodes
function get_node(mechanism::Mechanism{T,Nn,Ne,Nb}, id::Integer; origin::Bool=false) where {T,Nn,Ne,Nb}
    (origin && id == 0) && return mechanism.origin
    if id <= Ne
        return get_joint_constraint(mechanism, id)
    elseif id <= Ne+Nb
        return get_body(mechanism, id)
    else
        return get_contact_constraint(mechanism, id)
    end
end

function get_node(mechanism::Mechanism, name::Symbol)
    node = get_body(mechanism, name)
    if node === nothing
        node = get_joint_constraint(mechanism,name)
    end
    if node === nothing
        node = get_contact_constraint(mechanism,name)
    end
    return node
end

# maximal 
function get_maximal_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		x2 = body.state.x2
		v15 = body.state.v15
		q2 = body.state.q2
		ϕ15 = body.state.ϕ15
		pack_maximal_state!(z, x2, v15, q2, ϕ15, i)
	end
	return z
end

function get_maximal_gradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
    opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}

    step!(mechanism, z, u, opts=opts)
    jacobian_state, jacobian_control = get_maximal_gradients(mechanism)

    return jacobian_state, jacobian_control
end

# minimal 
function get_minimal_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	x = []

	mechanism = deepcopy(mechanism)
	timestep = mechanism.timestep

	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		c = zeros(T,0)
		v = zeros(T,0)
		pbody = get_body(mechanism, joint.parent_id)
		cbody = get_body(mechanism, joint.child_id)
		for (i, element) in enumerate([joint.translational, joint.rotational])
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

function get_minimal_velocities(mechanism::Mechanism{T}) where T
    d = Dict{Int,Vector{T}}()
    for joint in mechanism.joints
        push!(d, joint.id => minimal_velocities(mechanism, joint))
    end
    return d
end

function get_minimal_coordinates_velocities(mechanism::Mechanism{T}) where T
    d = Dict{Int,Vector{T}}()
    for joint in mechanism.joints
        push!(d, joint.id => [minimal_coordinates(mechanism, joint); minimal_velocities(mechanism, joint)])
    end
    return d
end

function get_minimal_configuration_vector(mechanism::Mechanism{T}) where T
    N = input_dimension(mechanism)
    x = zeros(T,N)
    off = 0
    for joint in mechanism.joints
        n = input_dimension(joint)
        x[off .+ (1:n)] += minimal_coordinates(mechanism, joint)
        off += n
    end
    return x
end

function get_minimal_velocity_vector(mechanism::Mechanism{T}) where T
    N = input_dimension(mechanism)
    x = zeros(T,N)
    off = 0
    for joint in mechanism.joints
        n = input_dimension(joint)
        x[off .+ (1:n)] += minimal_velocities(mechanism, joint)
        off += n
    end
    return x
end