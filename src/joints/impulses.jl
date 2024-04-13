function impulses!(mechanism, body::Body, joint::JointConstraint{T,Nλ}) where {T,Nλ}
    (Nλ > 0) && (body.state.d -= impulse_map(mechanism, joint, body) * joint.impulses[2])
    joint.spring && (body.state.d -= spring_impulses(mechanism, joint, body))
    joint.damper && (body.state.d -= damper_impulses(mechanism, joint, body))
    return
end

function impulse_map(mechanism, joint::JointConstraint, body::Body)
    relative = body.id == joint.parent_id ? :parent : :child
    pbody = get_body(mechanism, joint.parent_id)
    cbody = get_body(mechanism, joint.child_id)
    tra = impulse_map(relative, joint.translational,
        pbody, cbody,
        joint.impulses[2][joint_impulse_index(joint, 1)])
    rot = impulse_map(relative, joint.rotational,
        pbody, cbody,
        joint.impulses[2][joint_impulse_index(joint, 2)])
    return hcat(tra, rot)
end

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

################################################################################
# Impulse Transform
################################################################################
function impulse_transform(relative::Symbol, joint::Joint, xa::AbstractVector,
        qa::Quaternion, xb::AbstractVector, qb::Quaternion)
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=true)
    Diagonal([sones(3); 0.5 * sones(3)]) * transpose([X Q]) #TODO: 0.5 Q
end

################################################################################
# Derivatives
################################################################################
function impulse_map_jacobian(relative::Symbol, jacobian::Symbol, joint::Joint, pbody::Node{T}, cbody::Node{T}, λ) where T
    # ∂(G*λ)/∂(x,q)
    p = impulse_projector(joint) * λ
    impulse_transform_jacobian(relative, jacobian,
        joint,
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        p)
end

################################################################################
# Utilities
################################################################################
function get_joint_impulses(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    n1 = 1
    for j = 1:i-1
        n1 += impulses_length((joint.translational, joint.rotational)[j])
    end
    n2 = n1 - 1 + impulses_length((joint.translational, joint.rotational)[i])

    λi = SVector{n2-n1+1,T}(joint.impulses[2][SUnitRange(n1,n2)])
    return λi
end

function joint_impulse_index(joint::JointConstraint{T,N,Nc}, i::Int) where {T,N,Nc}
    s = 0
    for j = 1:i-1
        element = (joint.translational, joint.rotational)[j]
        s += impulses_length(element)
    end
    joint_impulse_index((joint.translational, joint.rotational)[i], s)
end