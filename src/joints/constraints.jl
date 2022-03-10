"""
    JointConstraint{T} <: Constraint{T}

    constraint restricting translational and rotational degrees of freedom between two Body objects.

    id: a unique identifying number
    name: a unique identifying name
    translational: Translational
    rotational: Rotational
    spring: flag for joint springs on
    damper: flag for joint dampers on
    parent_id: identifying number for parent Body{T}
    child_id: identifying number for child Body{T}
    minimal_index: indices for minimal coordinates
    impulses: joint impulses that maintain constraint between two Body{T} objects
"""
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

    function JointConstraint(data;
        name::Symbol=Symbol("joint_" * randstring(4)))

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
@generated function constraint(mechanism, joint::JointConstraint)
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    tra = :(constraint(joint.translational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint,1)], mechanism.μ, mechanism.timestep))
    rot = :(constraint(joint.rotational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint,2)], mechanism.μ, mechanism.timestep))
    return :(svcat($tra, $rot))
end

# constraints Jacobians
@generated function constraint_jacobian(joint::JointConstraint) 
    tra = :(constraint_jacobian(joint.translational, joint.impulses[2][joint_impulse_index(joint, 1)]))
    rot = :(constraint_jacobian(joint.rotational, joint.impulses[2][joint_impulse_index(joint, 2)]))
    return :(cat($tra, $rot, dims=(1,2)))
end

@generated function constraint_jacobian_configuration(mechanism, joint::JointConstraint, body::Body)
    relative = :(body.id == joint.parent_id ? :parent : :child)
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    tra = :(constraint_jacobian_configuration($relative,
        joint.translational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint, 1)], mechanism.timestep))
    rot = :(constraint_jacobian_configuration($relative,
        joint.rotational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint, 2)], mechanism.timestep))
    return :(vcat($tra, $rot))
end

# impulses
function impulses!(mechanism, body::Body, joint::JointConstraint)
    body.state.d -= impulse_map(mechanism, joint, body) * joint.impulses[2]
    joint.spring && (body.state.d -= spring_impulses(mechanism, joint, body))
    joint.damper && (body.state.d -= damper_impulses(mechanism, joint, body))
    return
end

@generated function impulse_map(mechanism, joint::JointConstraint, body::Body)
    relative = :(body.id == joint.parent_id ? :parent : :child)
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    tra = :(impulse_map($relative, joint.translational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint, 1)]))
    rot = :(impulse_map($relative, joint.rotational,
        $pbody, $cbody,
        joint.impulses[2][joint_impulse_index(joint, 2)]))
    return :(hcat($tra, $rot))
end

# impulses Jacobians
function impulses_jacobian_velocity!(mechanism, body::Body, joint::JointConstraint)

    # relative
    relative = (body.id == joint.parent_id ? :parent : (body.id == joint.child_id ? :child : error()))

    # time step
    timestep= mechanism.timestep

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
function off_diagonal_jacobians(mechanism, body::Body{T}, joint::JointConstraint)
    return -impulse_map(mechanism, joint, body), constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep)
end

function off_diagonal_jacobians(mechanism, joint::JointConstraint, body::Body)
    return constraint_jacobian_configuration(mechanism, joint, body) * integrator_jacobian_velocity(body, mechanism.timestep), -impulse_map(mechanism, joint, body)
end

function off_diagonal_jacobians(mechanism, pbody::Body, cbody::Body)
    timestep = mechanism.timestep

    dimpulse_map_parentb = szeros(6, 6)
    dimpulse_map_childa = szeros(6, 6)
    Ne = length(mechanism.joints)

    for connectionid in connections(mechanism.system, pbody.id)
        !(connectionid <= Ne) && continue # body
        joint = get_node(mechanism, connectionid)
        off = 0
        if pbody.id == joint.parent_id
            for element in (joint.translational, joint.rotational)
                Nj = length(element)
                if cbody.id == joint.child_id
                    joint.damper && (dimpulse_map_parentb -= damper_jacobian_velocity(:parent, :child, element, pbody, cbody, timestep))
                    joint.damper && (dimpulse_map_childa -= damper_jacobian_velocity(:child, :parent, element, pbody, cbody, timestep))
                end
                off += Nj
            end
        elseif cbody.id == joint.parent_id
            for element in (joint.translational, joint.rotational)
                Nj = length(element)
                if pbody.id == joint.child_id
                    joint.damper && (dimpulse_map_parentb -= damper_jacobian_velocity(:parent, :child, element, cbody, pbody, timestep))
                    joint.damper && (dimpulse_map_childa -= damper_jacobian_velocity(:child, :parent, element, cbody, pbody, timestep))
                end
                off += Nj
            end
        end
    end
    return dimpulse_map_parentb, dimpulse_map_childa
end

# linear system 
function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, joint::JointConstraint)
    matrix_entry.value = constraint_jacobian(joint)
    vector_entry.value = -constraint(mechanism, joint)
end

# springs
function spring_impulses(mechanism, joint::JointConstraint{T}, body::Body{T};
    unitary::Bool=false) where T

    relative = (body.id == joint.parent_id ? :parent : :child)
    impulses = szeros(T,6)

    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)

    impulses += spring_impulses(relative, joint.translational,
        pbody,
        cbody, 
        mechanism.timestep,
        unitary=unitary)

    impulses += spring_impulses(relative, joint.rotational,
        pbody, 
        cbody,
        mechanism.timestep,
        unitary=unitary)

    return impulses
end

# dampers
function damper_impulses(mechanism, joint::JointConstraint{T}, body::Body;
    unitary::Bool=false) where T

    relative = (body.id == joint.parent_id ? :parent : :child)
    impulses = szeros(T,6)

    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)

    impulses += damper_impulses(relative, joint.translational,
        pbody, 
        cbody,
        mechanism.timestep,
        unitary=unitary)

    impulses += damper_impulses(relative, joint.rotational,
        pbody,
        cbody,
        mechanism.timestep,
        unitary=unitary)
        
    return impulses
end

# inputs
function set_input!(joint::JointConstraint{T,N,Nc}, input::AbstractVector) where {T,N,Nc}
    @assert length(input) == input_dimension(joint)
    # translational
    r_idx = SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])
    length(r_idx) > 0 && set_input!(joint.translational, input[r_idx])
    # rotational
    r_idx = SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])
    length(r_idx) > 0 && set_input!(joint.rotational, input[r_idx])
    return
end

function add_input!(joint::JointConstraint{T,N,Nc}, input::AbstractVector) where {T,N,Nc}
    @assert length(input) == input_dimension(joint)
    add_input!(joint.translational, input[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])])
    add_input!(joint.rotational, input[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])])
    return
end

@generated function input_jacobian_control(mechanism, joint::JointConstraint{T,N,Nc}, body::Body) where {T,N,Nc}
    relative = :(body.id == joint.parent_id ? :parent : :child)
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    rot = :(input_jacobian_control($relative, joint.translational, $pbody, $cbody))
    tra = :(input_jacobian_control($relative, joint.rotational, $pbody, $cbody))
    return :(hcat($rot, $tra))
end

function input_impulse!(joint::JointConstraint{T,N,Nc}, mechanism, clear::Bool=true) where {T,N,Nc}
    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)
    input_impulse!(joint.translational, pbody, cbody, mechanism.timestep, clear)
    input_impulse!(joint.rotational, pbody, cbody, mechanism.timestep, clear)
    return
end

# minimal
@generated function minimal_coordinates(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    tra = :(minimal_coordinates(joint.translational, $pbody, $cbody))
    rot = :(minimal_coordinates(joint.rotational, $pbody, $cbody))
    return :(svcat($tra, $rot))
end

@generated function minimal_velocities(mechanism, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    pbody = :(get_body(mechanism, joint.parent_id))
    cbody = :(get_body(mechanism, joint.child_id))
    tra = :(minimal_velocities(joint.translational, $pbody, $cbody, mechanism.timestep))
    rot = :(minimal_velocities(joint.rotational, $pbody, $cbody, mechanism.timestep))
    return :(svcat($tra, $rot))
end

################################################################################
# Utilities
################################################################################
function get_joint_impulses(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    n1 = 1
    for j = 1:i-1
        n1 += impulses_length([joint.translational, joint.rotational][j])
    end
    n2 = n1 - 1 + impulses_length([joint.translational, joint.rotational][i])

    λi = SVector{n2-n1+1,T}(joint.impulses[2][n1:n2])
    return λi
end

function joint_impulse_index(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    s = 0
    for j = 1:i-1
        element = [joint.translational, joint.rotational][j]
        s += impulses_length(element)
    end
    joint_impulse_index([joint.translational, joint.rotational][i], s) # to be allocation free
end

function reset!(joint::JointConstraint{T,N,Nc};
    scale::T=1.0) where {T,N,Nc}
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

function set_spring_damper_values!(joints, spring, damper;
    ignore_origin::Bool=true)
    i = 1
    for joint in joints
        (ignore_origin && joint.parent_id == 0) && continue
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

function input_dimension(joint::JointConstraint{T,N,Nc};
    ignore_floating_base::Bool=false) where {T,N,Nc}
    ignore_floating_base && (N == 0) && return 0
    N̄ = 0
    N̄ = input_dimension(joint.translational) + input_dimension(joint.rotational)
    return N̄
end
