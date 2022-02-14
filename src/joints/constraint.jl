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
        T = getT(data[1][1])

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
            Nλ = λlength(joint)
            Nset = ηlength(joint)
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

# TODO currently assumed constraints are in order and only joints which is the case unless very low level constraint setting
function set_joint_position!(mechanism, joint::JointConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    Nλ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ += λlength(element)
    end
    @assert length(xθ)==3*Nc-Nλ

    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    Δx = xθ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δθ = xθ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    set_minimal_coordinates!(body1, body2, joint, Δx=Δx, Δθ=Δθ)
    return body2.state.x2[1], body2.state.q2[1]
end

function set_velocity!(mechanism, joint::JointConstraint{T,N,Nc}, vϕ) where {T,N,Nc}
    Nλ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ += λlength(element)
    end

    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    Δv = vϕ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δϕ = vϕ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    set_minimal_velocities_new!(body1, body2, joint, mechanism.timestep, Δv=Δv, Δϕ=Δϕ)
    return body2.state.v15, body2.state.ϕ15
end

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

"""
    minimal_coordinates(mechanism, jointonstraint)

Gets the minimal coordinates of joint `jointonstraint`.
"""
@generated function minimal_coordinates(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimal_coordinates([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id))) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@generated function minimal_velocities_new(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(minimal_velocities_new([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), mechanism.timestep)) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@inline function impulses!(mechanism, body::Body, joint::JointConstraint)
    body.state.d -= impulse_map(mechanism, joint, body) * joint.impulses[2]
    joint.spring && (body.state.d -= apply_spring(mechanism, joint, body))
    joint.damper && (body.state.d -= apply_damper(mechanism, joint, body))
    return
end

@inline function impulses_jacobian_velocity!(mechanism, body::Body, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    # @warn "maybe need some work"
    if body.id == joint.parent_id
        impulses_jacobian_parent!(mechanism, body, joint)
    elseif body.id == joint.child_id
        impulses_jacobian_child!(mechanism, body, joint)
    else
        error()
    end
    return
end

function impulses_jacobian_parent!(mechanism, pbody::Body, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    timestep = mechanism.timestep
    _, _, q2, ω2 = current_configuration_velocity(pbody.state)
    M = integrator_jacobian_velocity(q2, ω2, timestep)

    # child boy
    cbody = get_body(mechanism, joint.child_id)

    off = 0
    for element in [joint.translational, joint.rotational]
        Nj = length(element)
        joint.spring && (pbody.state.D -= spring_parent_jacobian_velocity_parent(element, pbody, cbody, timestep))
        joint.damper && (pbody.state.D -= damper_parent_jacobian_velocity_parent(element, pbody, cbody, timestep))
        off += Nj
    end
    return nothing
end

function impulses_jacobian_child!(mechanism, cbody::Body, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    timestep = mechanism.timestep
    x2, v2, q2, ω2 = current_configuration_velocity(cbody.state)
    M = integrator_jacobian_velocity(q2, ω2, timestep)

    off = 0
    for element in [joint.translational, joint.rotational]
        if joint.child_id == cbody.id
            pbody = get_body(mechanism, joint.parent_id)
            joint.spring && (cbody.state.D -= spring_child_jacobian_velocity_child(element, pbody, cbody, timestep))
            joint.damper && (cbody.state.D -= damper_child_jacobian_velocity_child(element, pbody, cbody, timestep))
        end
    end
    return nothing
end

@generated function constraint(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(constraint([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), joint.impulses[2][λindex(joint,$i)], mechanism.timestep)) for i = 1:Nc]
    return :(svcat($(vec...)))
end

@inline function apply_spring(mechanism, joint::JointConstraint, body::Body; unitary::Bool=false)
    body.id == joint.parent_id ? (return spring_parent(mechanism, joint, body, unitary=unitary)) : (return spring_child(mechanism, joint, body, unitary=unitary))
end

@inline function apply_damper(mechanism, joint::JointConstraint, body::Body; unitary::Bool=false)
    body.id == joint.parent_id ? (return damper_parent(mechanism, joint, body, unitary=unitary)) : (return damper_child(mechanism, joint, body, unitary=unitary))
end

@inline function off_diagonal_jacobians(mechanism, body::Body{T}, joint::JointConstraint{T,N}) where {T,N}
    return -impulse_map(mechanism, joint, body), constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep)
end

@inline function off_diagonal_jacobians(mechanism, joint::JointConstraint{T,N}, body::Body{T}) where {T,N}
    return constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep), -impulse_map(mechanism, joint, body)
end

@generated function impulse_map_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(impulse_map_parent([joint.translational, joint.rotational][$i], body, get_body(mechanism, joint.child_id), joint.child_id, joint.impulses[2][λindex(joint,$i)], mechanism.timestep)) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@generated function impulse_map_child(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(impulse_map_child([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), body, joint.child_id, joint.impulses[2][λindex(joint,$i)], mechanism.timestep)) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@generated function constraint_jacobian_configuration(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    vec = [:(constraint_jacobian_configuration([joint.translational, joint.rotational][$i], joint.impulses[2][λindex(joint,$i)])) for i = 1:Nc]
    return :(cat($(vec...), dims=(1,2)))
end

@generated function constraint_jacobian_configuration(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    # relatives = (body.id == joint.parent_id ? :parent : :child)
    vec = [:(constraint_jacobian_configuration(body.id == joint.parent_id ? :parent : :child, [joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), joint.child_id, joint.impulses[2][λindex(joint,$i)], mechanism.timestep)) for i = 1:Nc]
    return :(vcat($(vec...)))
end

@inline function spring_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += spring_parent([joint.translational, joint.rotational][i], body, get_body(mechanism, joint.child_id), mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

@inline function spring_child(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += spring_child([joint.translational, joint.rotational][i], get_body(mechanism, joint.parent_id), body, mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

@inline function damper_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damper_parent([joint.translational, joint.rotational][i], body, get_body(mechanism, joint.child_id), mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

@inline function damper_child(mechanism, joint::JointConstraint{T,N,Nc}, body::Body; unitary::Bool=false) where {T,N,Nc}
    vec = szeros(T,6)
    for i=1:Nc
        vec += damper_child([joint.translational, joint.rotational][i], get_body(mechanism, joint.parent_id), body, mechanism.timestep, joint.child_id, unitary=unitary)
    end
    return vec
end

@inline function input_jacobian_control(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    body.id == joint.parent_id ? (return input_jacobian_control_parent(mechanism, joint, body)) : (return input_jacobian_control_child(mechanism, joint, body))
end

@generated function input_jacobian_control_parent(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(input_jacobian_control_parent([joint.translational, joint.rotational][$i], body, get_body(mechanism, joint.child_id), joint.child_id)) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@generated function input_jacobian_control_child(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    vec = [:(input_jacobian_control_child([joint.translational, joint.rotational][$i], get_body(mechanism, joint.parent_id), body, joint.child_id)) for i = 1:Nc]
    return :(hcat($(vec...)))
end

@inline function apply_input!(joint::JointConstraint{T,N,Nc}, mechanism, clear::Bool=true) where {T,N,Nc}
    for i=1:Nc
        apply_input!([joint.translational, joint.rotational][i], get_body(mechanism, joint.parent_id), get_body(mechanism, joint.child_id), mechanism.timestep, clear)
    end
    return
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

@inline function impulse_map(mechanism, constraint::Constraint, body::Body)
    if body.id == constraint.parent_id
        return impulse_map_parent(mechanism, constraint, body)
    else
        return impulse_map_child(mechanism, constraint, body)
    end
end

function λindex(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    s = 0
    for j = 1:i-1
        element = [joint.translational, joint.rotational][j]
        s += ηlength(element)
    end
    λindex([joint.translational, joint.rotational][i], s) # to be allocation free
end

function reset!(joint::JointConstraint{T,N,Nc}; scale::T=1.0) where {T,N,Nc}
    λ = []
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ = λlength(element)
        Nb = blength(element)
        push!(λ, [scale * sones(2Nb); szeros(Nλ)])
    end
    joint.impulses[1] = vcat(λ...)
    joint.impulses[2] = vcat(λ...)
    return
end
