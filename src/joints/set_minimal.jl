
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
        Δx::AbstractVector=szeros(control_dimension(joint.constraints[1])),
        Δθ::AbstractVector=szeros(control_dimension(joint.constraints[2])),
        Δv::AbstractVector=szeros(control_dimension(joint.constraints[1])),
        Δϕ::AbstractVector=szeros(control_dimension(joint.constraints[2])))
    # We need to set the minimal coordinates of the rotational joint first
    # since xb = fct(qb, Δx)
    # since vb = fct(ϕb, Δv)
    set_minimal_coordinates!(pnode, cnode, joint; Δx=Δx, Δθ=Δθ)
    set_minimal_velocities!(pnode, cnode, joint; Δv=Δv, Δϕ=Δϕ)
    return nothing
end

function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::JointConstraint;
        Δx::AbstractVector=szeros(control_dimension(joint.constraints[1])),
        Δθ::AbstractVector=szeros(control_dimension(joint.constraints[2])))
    # We need to set the minimal coordinates of the rotational joint first
    # since xb = fct(qb, Δx)
    set_minimal_coordinates!(pnode, cnode, joint.constraints[2]; Δθ=Δθ)
    set_minimal_coordinates!(pnode, cnode, joint.constraints[1]; Δx=Δx)
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
function set_minimal_velocities!(pnode::Node, cnode::Node, joint::JointConstraint;
        Δv::AbstractVector=szeros(control_dimension(joint.constraints[1])),
        Δϕ::AbstractVector=szeros(control_dimension(joint.constraints[2])))
    # We need to set the minimal coordinates of the rotational joint first
    # since vb = fct(ϕb, Δv)
    set_minimal_velocities!(pnode, cnode, joint.constraints[2]; Δϕ=Δϕ)
    set_minimal_velocities!(pnode, cnode, joint.constraints[1]; Δv=Δv)
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
