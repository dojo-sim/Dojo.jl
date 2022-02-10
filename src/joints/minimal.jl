################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::Joint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::JointConstraint;
        Δx::AbstractVector=szeros(control_dimension(joint.translational)),
        Δθ::AbstractVector=szeros(control_dimension(joint.rotational)))
    # We need to set the minimal coordinates of the rotational joint first
    # since xb = fct(qb, Δx)
    set_minimal_coordinates!(pnode, cnode, joint.rotational; Δθ=Δθ)
    set_minimal_coordinates!(pnode, cnode, joint.translational; Δx=Δx)
    return nothing
end

function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::Rotational;
        Δθ::AbstractVector=szeros(control_dimension(joint)))
        # Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
    qoffset = joint.qoffset
    qa = pnode.state.q2[1]
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δq = axis_angle_to_quaternion(Aᵀ*Δθ)
    qb = qa * qoffset * Δq
    set_position!(cnode; x=cnode.state.x2[1], q = qb)
    return nothing
end

function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::Translational;
        Δx::AbstractVector=szeros(control_dimension(joint)))
        # Δx is expressed in along the joint's nullspace axes, in pnode's frame

    pa = joint.vertices[1]
    pb = joint.vertices[2]

    qa = pnode.state.q2[1]
    xa = pnode.state.x2[1]

    qb = cnode.state.q2[1]

    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    xb = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
    set_position!(cnode; x = xb, q=cnode.state.q2[1])
    return nothing
end

################################################################################
# Velocities
################################################################################
@inline function minimal_velocities(joint::Joint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.x2[1], statea.v15, statea.q2[1], statea.ϕ15,
		stateb.x2[1], stateb.v15, stateb.q2[1], stateb.ϕ15)
end

function minimal_velocities_jacobian_configuration(jacobian_relative::Symbol,
        joint::Joint, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ϕa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    (jacobian_relative == :parent) && (return minimal_velocities_jacobian_configuration_parent(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    (jacobian_relative == :child) && (return minimal_velocities_jacobian_configuration_child(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    return
end

function minimal_velocities_jacobian_velocity(jacobian_relative::Symbol,
        joint::Joint, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ϕa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    (jacobian_relative == :parent) && (return minimal_velocities_jacobian_velocity_parent(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    (jacobian_relative == :child) && (return minimal_velocities_jacobian_velocity_child(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    return
end

function set_minimal_velocities!(pnode::Node, cnode::Node, joint::JointConstraint;
        Δv::AbstractVector=szeros(control_dimension(joint.translational)),
        Δϕ::AbstractVector=szeros(control_dimension(joint.rotational)))
    # We need to set the minimal coordinates of the rotational joint first
    # since vb = fct(ϕb, Δv)
    set_minimal_velocities!(pnode, cnode, joint.rotational; Δϕ=Δϕ)
    set_minimal_velocities!(pnode, cnode, joint.translational; Δv=Δv)
    return nothing
end

function set_minimal_velocities!(pnode::Node, cnode::Node, joint::Rotational;
        Δϕ::AbstractVector=szeros(control_dimension(joint)))
        # Δϕ is expressed in along the joint's nullspace axes, in pnode's offset frame
        # We need to set the minimal coordinates of the rotational joint first
        # since ϕb = fct(qb, Δϕ)
    qoffset = joint.qoffset
    qa = pnode.state.q2[1]
    qb = cnode.state.q2[1]
    ϕa = pnode.state.ϕ15
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δϕ_a = vrotate(Aᵀ*Δϕ, qoffset) # in pnode's frame
    ϕb = vrotate(ϕa + Δϕ_a , inv(qb) * qa)
    set_velocity!(cnode; v=cnode.state.v15, ω=ϕb)
    return nothing
end

function set_minimal_velocities!(pnode::Node, cnode::Node, joint::Translational;
        Δv::AbstractVector=szeros(control_dimension(joint)))
        # Δv is expressed in along the joint's nullspace axes, in pnode's frame
    xa = pnode.state.x2[1]
    va = pnode.state.v15
    qa = pnode.state.q2[1]
    ϕa = pnode.state.ϕ15

    xb = cnode.state.x2[1]
    vb = cnode.state.v15
    qb = cnode.state.q2[1]
    ϕb = cnode.state.ϕ15

    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δv = Aᵀ * Δv # in pnode's frame
    Δvw = vrotate(Δv, qa) # in the world frame
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    vb = Δvw - skew(pbcb_w) * vrotate(ϕb, qb) + (va + skew(pbca_w) * vrotate(ϕa, qa)) # in world frame
    set_velocity!(cnode; v=vb, ω=cnode.state.ϕ15)
    return nothing
end

################################################################################
# Coordinates and Velocities
################################################################################
function set_minimal_coordinates_velocities!(mechanism::Mechanism, joint::JointConstraint;
        xmin::AbstractVector=szeros(2control_dimension(joint)))
    pnode = get_body(mechanism, joint.parent_id)
    cnode = get_body(mechanism, joint.child_id)
    nu = control_dimension(joint)
    Δx = xmin[SUnitRange(joint.minimal_index[1]...)]
    Δθ = xmin[SUnitRange(joint.minimal_index[2]...)]
    Δv = xmin[nu .+ SUnitRange(joint.minimal_index[1]...)]
    Δϕ = xmin[nu .+ SUnitRange(joint.minimal_index[2]...)]
    set_minimal_coordinates_velocities!(pnode, cnode, joint; Δx=Δx, Δθ=Δθ, Δv=Δv, Δϕ=Δϕ)
end

function set_minimal_coordinates_velocities!(pnode::Node, cnode::Node, joint::JointConstraint;
        Δx::AbstractVector=szeros(control_dimension(joint.translational)),
        Δθ::AbstractVector=szeros(control_dimension(joint.rotational)),
        Δv::AbstractVector=szeros(control_dimension(joint.translational)),
        Δϕ::AbstractVector=szeros(control_dimension(joint.rotational)))
    # We need to set the minimal coordinates of the rotational joint first
    # since xb = fct(qb, Δx)
    # since vb = fct(ϕb, Δv)
    set_minimal_coordinates!(pnode, cnode, joint; Δx=Δx, Δθ=Δθ)
    set_minimal_velocities!(pnode, cnode, joint; Δv=Δv, Δϕ=Δϕ)
    return nothing
end
