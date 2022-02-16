mutable struct JointConstraint{T,N,Nc,TJ,RJ} <: Constraint{T,N}
    # ID
    id::Int64
    name::Symbol

    # joint constraints
    translational::TJ
    rotational::RJ

    # springs and dampers
    spring::Bool
    damper::Bool

    # neighbor IDs
    parent_id::Int
    child_id::Int

    # indices
    minimal_index::SVector{Nc,SVector{2,Int64}} # indices for minimal coordinates, assumes joints # Nc = 2 THIS IS SPECIAL CASED

    # impulses
    impulses::Vector{SVector{N,T}}

    function JointConstraint(data; name::Symbol=Symbol("joint_" * randstring(4)))
        @assert data[1][2] == data[2][2] # check parent ids
        @assert data[1][3] == data[2][3] # check child ids

        # joints
        translational = data[1][1]
        rotational = data[2][1]

        # IDs
        parent_id = data[1][2]
        child_id = data[1][3]

        # data dype
        T = typeof(data[1][1]).parameters[1]

        # set springs & dampers off
        spring = false
        damper = false

        minimal_index = Vector{Int64}[]
        N = 0
        for joint_data in data
            joint = joint_data[1]

            # set spring & damper on
            joint.spring != 0 && (spring = true)
            joint.damper != 0 && (damper = true)

            # minimal-coordaintes indices
            Nλ = joint_length(joint)
            Nset = impulses_length(joint)
            if isempty(minimal_index)
                push!(minimal_index, [1;3-Nλ])
            else
                push!(minimal_index, [last(minimal_index)[2]+1; last(minimal_index)[2]+3-Nλ])
            end
            N += Nset
        end

        Nc = 2
        impulses = [zeros(T, N) for i=1:2]

        return new{T,N,Nc,typeof(translational),typeof(rotational)}(getGlobalID(), name, translational, rotational, spring, damper, parent_id, child_id, minimal_index, impulses)
    end
end

# constraints
@generated function constraint(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    tra = :(constraint(joint.translational, 
        get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id),
        joint.impulses[2][joint_impulse_index(joint,1)], mechanism.timestep))
    rot = :(constraint(joint.rotational,    
        get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), 
        joint.impulses[2][joint_impulse_index(joint,2)], mechanism.timestep))
    return :(svcat($tra, $rot))
end

# constraints Jacobians 
@generated function constraint_jacobian_configuration(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    tra = :(constraint_jacobian_configuration(joint.translational, joint.impulses[2][joint_impulse_index(joint, 1)]))
    rot = :(constraint_jacobian_configuration(joint.rotational, joint.impulses[2][joint_impulse_index(joint, 2)]))
    return :(cat($tra, $rot, dims=(1,2)))
end

@generated function constraint_jacobian_configuration(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    relative = :(body.id == joint.parent_id ? :parent : :child)
    tra = :(constraint_jacobian_configuration($relative, 
        joint.translational, 
        get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), 
        joint.child_id, joint.impulses[2][joint_impulse_index(joint, 1)], mechanism.timestep))
    rot = :(constraint_jacobian_configuration($relative, 
        joint.rotational, 
        get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), 
        joint.child_id, joint.impulses[2][joint_impulse_index(joint, 2)], mechanism.timestep))
    return :(vcat($tra, $rot))
end

# impulses 
@inline function impulses!(mechanism, body::Body, joint::JointConstraint)
    body.state.d -= impulse_map(mechanism, joint, body) * joint.impulses[2]
    joint.spring && (body.state.d -= apply_spring(mechanism, joint, body))
    joint.damper && (body.state.d -= apply_damper(mechanism, joint, body))
    return
end

@generated function impulse_map(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    relative = :(body.id == joint.parent_id ? :parent : :child)
    tra = :(impulse_map($relative, joint.translational, 
        get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), 
        joint.child_id, joint.impulses[2][joint_impulse_index(joint, 1)], mechanism.timestep))
    rot = :(impulse_map($relative, joint.rotational,
        get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), 
        joint.child_id, joint.impulses[2][joint_impulse_index(joint, 2)], mechanism.timestep))
    return :(hcat($tra, $rot))
end

# impulses Jacobians
@inline function impulses_jacobian_velocity!(mechanism, body::Body, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}

    # relative
    relative = (body.id == joint.parent_id ? :parent : (body.id == joint.child_id ? :child : error()))

    # time step
    timestep = mechanism.timestep
    
    # bodies
    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)

    # springs
    joint.spring && (body.state.D -= spring_jacobian_velocity(relative, relative, joint.translational, pbody, cbody, timestep))
    joint.spring && (body.state.D -= spring_jacobian_velocity(relative, relative, joint.rotational, pbody, cbody, timestep))

    # dampers
    joint.damper && (body.state.D -= damper_jacobian_velocity(relative, relative, joint.translational, pbody, cbody, timestep))
    joint.damper && (body.state.D -= damper_jacobian_velocity(relative, relative, joint.rotational, pbody, cbody, timestep))

    return
end

# off-diagonal Jacobians 
@inline function off_diagonal_jacobians(mechanism, body::Body{T}, joint::JointConstraint{T,N}) where {T,N}
    return -impulse_map(mechanism, joint, body), constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep)
end

@inline function off_diagonal_jacobians(mechanism, joint::JointConstraint{T,N}, body::Body{T}) where {T,N}
    return constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep), -impulse_map(mechanism, joint, body)
end

# springs 
@inline function apply_spring(mechanism, joint::JointConstraint, body::Body; unitary::Bool=false)
    body.id == joint.parent_id ? (return spring_parent(mechanism, joint, body, unitary=unitary)) : (return spring_child(mechanism, joint, body, unitary=unitary))
end

@inline function spring_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += spring_force(:parent, [joint.translational, joint.rotational][i], body, get_body(mechanism, joint.child_id), mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

@inline function spring_child(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += spring_force(:child, [joint.translational, joint.rotational][i], get_body(mechanism, joint.parent_id), body, mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

# dampers
@inline function apply_damper(mechanism, joint::JointConstraint, body::Body; unitary::Bool=false)
    body.id == joint.parent_id ? (return damper_parent(mechanism, joint, body, unitary=unitary)) : (return damper_child(mechanism, joint, body, unitary=unitary))
end

@inline function damper_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damper_force(:parent, [joint.translational, joint.rotational][i], body, get_body(mechanism, joint.child_id), mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

@inline function damper_child(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damper_force(:child, [joint.translational, joint.rotational][i], get_body(mechanism, joint.parent_id), body, mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

# inputs 
function set_input!(joint::JointConstraint{T,N,Nc}, input::AbstractVector) where {T,N,Nc}
    @assert length(input)==control_dimension(joint)
    for i = 1:Nc
        r_idx = SUnitRange(joint.minimal_index[i][1], joint.minimal_index[i][2])
        length(r_idx) == 0 && continue
        set_input!([joint.translational, joint.rotational][i], input[SUnitRange(joint.minimal_index[i][1], joint.minimal_index[i][2])])
    end
    return
end

function add_input!(joint::JointConstraint{T,N,Nc}, input::AbstractVector) where {T,N,Nc}
    @assert length(input)==control_dimension(joint)
    for i = 1:Nc
        add_input!([joint.translational, joint.rotational][i], input[SUnitRange(joint.minimal_index[i][1], joint.minimal_index[i][2])])
    end
    return
end

@generated function input_jacobian_control(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(input_jacobian_control((body.id == joint.parent_id ? :parent : :child), [joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), joint.child_id)) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@inline function input_impulse!(joint::JointConstraint{T,N,Nc}, mechanism, clear::Bool=true) where {T,N,Nc}
    for i=1:Nc
        input_impulse!([joint.translational, joint.rotational][i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), mechanism.timestep, clear)
    end
    return
end

# set methods
function set_position!(mechanism, joint::JointConstraint, xθ; iter::Bool=true)
    if !iter
        set_joint_position!(mechanism, joint, xθ)
    else
        currentvals = minimal_coordinates(mechanism)
        set_joint_position!(mechanism, joint, xθ)
        for id in recursivedirectchildren!(mechanism.system, joint.id)
            node = get_node(mechanism, id)
            if node isa JointConstraint
                set_joint_position!(mechanism, node, currentvals[id])
            end
        end
    end
    return
end

function set_joint_position!(mechanism, joint::JointConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    Nλ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ += joint_length(element)
    end
    @assert length(xθ)==3*Nc-Nλ

    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    Δx = xθ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δθ = xθ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    set_minimal_coordinates!(body1, body2, joint, mechanism.timestep, Δx=Δx, Δθ=Δθ)
    return body2.state.x2[1], body2.state.q2[1]
end

function set_velocity!(mechanism, joint::JointConstraint{T,N,Nc}, vϕ) where {T,N,Nc}
    Nλ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ += joint_length(element)
    end

    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    Δv = vϕ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δϕ = vϕ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    set_minimal_velocities!(body1, body2, joint, mechanism.timestep, Δv=Δv, Δϕ=Δϕ)
    return body2.state.v15, body2.state.ϕ15
end

# minimal
@generated function minimal_coordinates(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimal_coordinates([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@generated function minimal_velocities(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimal_velocities([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), mechanism.timestep)) for i = 1:Nc]
    return :(svcat($(vec...)))
end

function set_spring_damper_values!(joints, spring, damper)
    i = 1
    for joint in joints
        joint.parent_id == 0 && continue
        k = (length(spring) > 1) ? spring[i] : spring
        b = (length(damper) > 1) ? damper[i] : damper
        joint.spring = k > 0.0
        joint.damper = b > 0.0
        for element in [joint.translational, joint.rotational]
            element.spring = max(0.0, k)
            element.damper = max(0.0, b)
        end
        i += 1
    end
    return joints
end

################################################################################
# Utilities
################################################################################
function joint_impulse_index(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    s = 0
    for j = 1:i-1
        element = [joint.translational, joint.rotational][j]
        s += impulses_length(element)
    end
    joint_impulse_index([joint.translational, joint.rotational][i], s) # to be allocation free
end

function reset!(joint::JointConstraint{T,N,Nc}; scale::T=1.0) where {T,N,Nc}
    λ = []
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ = joint_length(element)
        Nb = limits_length(element)
        push!(λ, [scale * sones(2Nb); szeros(Nλ)])
    end
    joint.impulses[1] = vcat(λ...)
    joint.impulses[2] = vcat(λ...)
    return
end
