################################################################################
# SET: Coordinates Joint constraints
################################################################################
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

function set_minimal_coordinates_jacobian_parent(pnode::Node, cnode::Node, joint::Rotational{T};
    Δθ::AbstractVector=szeros(control_dimension(joint))) where T
    # Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
    # qoffset = joint.qoffset
    # qa = pnode.state.q2[1]
    # Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    # Δq = axis_angle_to_quaternion(Aᵀ*Δθ)
    # qb = qa * qoffset * Δq
    # set_position!(cnode; x=cnode.state.x2[1], q = qb)
    # [
    #     cnode.state.x2[1] = self;
    #     cnode.state.q2[1] = qa * qoffset * Δq
    # ]
    X = SMatrix{3,3,T,9}(Diagonal(szeros(3)))
    Q = Rmat(qoffset * Δq)
    return X, Q
end

function set_minimal_coordinates_jacobian_minimal(pnode::Node, cnode::Node, joint::Rotational{T};
    Δθ::AbstractVector=szeros(control_dimension(joint))) where T
    # Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
    # qoffset = joint.qoffset
    # qa = pnode.state.q2[1]
    # Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    # Δq = axis_angle_to_quaternion(Aᵀ*Δθ)
    # qb = qa * qoffset * Δq
    # set_position!(cnode; x=cnode.state.x2[1], q = qb)
    # [
    #     x = self;
    #     q = qa * qoffset * Δq
    # ]

    nu = control_dimension(joint)
    xθ = SMatrix{3,nu,T,3 * nu}(zeros(3, nu))

    ∂axis_angle_to_quaternion∂a = FiniteDiff.finite_difference_jacobian(r -> axis_angle_to_quaternion(r), Aᵀ * Δθ)
    qθ = Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂a * Aᵀ

    return [xθ; qθ]
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
    xb = xa + vrotate(pa + Aᵀ * Δx, qa) - vrotate(pb, qb)
    set_position!(cnode; x = xb, q=cnode.state.q2[1])
    return nothing
end

function set_minimal_coordinates_jacobian_parent(pnode::Node, cnode::Node, joint::Translational;
    Δx::AbstractVector=szeros(control_dimension(joint)))
    # Δx is expressed in along the joint's nullspace axes, in pnode's frame

    # pa = joint.vertices[1]
    # pb = joint.vertices[2]

    # qa = pnode.state.q2[1]
    # xa = pnode.state.x2[1]

    # qb = cnode.state.q2[1]

    # Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    # xb = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
    # set_position!(cnode; x = xb, q=cnode.state.q2[1])
    # [
    #     x = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
    #     q = self
    # ]
    xx = SMatrix{3,3,T,9}(Diagonal(sones(3)))
    xq = ∂vrotate∂q(pa + Aᵀ*Δx, qa)
    qx = SMatrix{4,3,T,12}(Diagonal(szeros(3)))
    qq = SMatrix{4,4,T,16}(Diagonal(szeros(4)))

    return [xx xq; qx qq]
end

function set_minimal_coordinates_jacobian_minimal(pnode::Node, cnode::Node, joint::Translational;
    Δx::AbstractVector=szeros(control_dimension(joint)))
    # Δx is expressed in along the joint's nullspace axes, in pnode's frame

    # pa = joint.vertices[1]
    # pb = joint.vertices[2]

    # qa = pnode.state.q2[1]
    # xa = pnode.state.x2[1]

    # qb = cnode.state.q2[1]

    # Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    # xb = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
    # set_position!(cnode; x = xb, q=cnode.state.q2[1])
    # [
    #     x = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
    #     q = self
    # ]

    nu = control_dimension(joint)

    xθ = ∂vrotate∂p(pa + Aᵀ * Δx) * Aᵀ
    qθ = SMatrix{4,nu,T,4 * nu}(zeros(4, nu))

    return [xθ; qθ]
end

################################################################################
# SET: Velocities Joint constraints
################################################################################
function set_minimal_velocities_new!(pnode::Node, cnode::Node, joint::JointConstraint, timestep;
		Δx::AbstractVector=szeros(control_dimension(joint.translational)),
		Δθ::AbstractVector=szeros(control_dimension(joint.rotational)),
        Δv::AbstractVector=szeros(control_dimension(joint.translational)),
        Δϕ::AbstractVector=szeros(control_dimension(joint.rotational)))
    # We need to set the minimal coordinates of the rotational joint first
    # since vb = fct(ϕb, Δv)
    # set_minimal_velocities!(pnode, cnode, joint.rotational; Δϕ=Δϕ)
    # set_minimal_velocities!(pnode, cnode, joint.translational; Δv=Δv)
	rot = joint.rotational
	tra = joint.translational
	pa = tra.vertices[1]
    pb = tra.vertices[2]
	qoffset = rot.qoffset
	Arotᵀ = zerodimstaticadjoint(nullspace_mask(rot))
	Atraᵀ = zerodimstaticadjoint(nullspace_mask(tra))

	xa, va, qa, ϕa = initial_configuration_velocity(pnode.state)
	xb, qb = current_configuration(cnode.state)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)

	# Finite difference
	Δx1 = Δx - Δv * timestep
	Δθ1 = Δθ - Δϕ * timestep

	# Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
    Δq1 = axis_angle_to_quaternion(Arotᵀ*Δθ1)
    qb1 = qa1 * qoffset * Δq1
    # Δx is expressed in along the joint's nullspace axes, in pnode's frame
    xb1 = xa1 + vrotate(pa + Atraᵀ*Δx1, qa1) - vrotate(pb, qb1)

	# Finite difference
	vb = (xb - xb1) / timestep
	ϕb = angular_velocity(qb1, qb, timestep)

	set_velocity!(cnode; v=vb, ω=ϕb)
    return nothing
end


# end

# function set_minimal_coordinates_jacobian_minimal(pnode::Node, cnode::Node, joint::Translational;
#     Δx::AbstractVector=szeros(control_dimension(joint)))
#     # Δx is expressed in along the joint's nullspace axes, in pnode's frame
#
#     # pa = joint.vertices[1]
#     # pb = joint.vertices[2]
#
#     # qa = pnode.state.q2[1]
#     # xa = pnode.state.x2[1]
#
#     # qb = cnode.state.q2[1]
#
#     # Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
#     # xb = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
#     # set_position!(cnode; x = xb, q=cnode.state.q2[1])
#     # [
#     #     x = xa + vrotate(pa + Aᵀ*Δx, qa) - vrotate(pb, qb)
#     #     q = self
#     # ]
#
#     nu = control_dimension(joint)
#
#     xθ = ∂vrotate∂p(pa + Aᵀ * Δx) * Aᵀ
#     qθ = SMatrix{4,nu,T,4 * nu}(zeros(4, nu))
#
#     return [xθ; qθ]
# end
#
# ################################################################################
# # Velocities
# ################################################################################
# @inline function minimal_velocities(joint::Joint, body1::Node, body2::Node)
#     statea = body1.state
#     stateb = body2.state
#     return minimal_velocities(joint, statea.x2[1], statea.v15, statea.q2[1], statea.ϕ15,
# 		stateb.x2[1], stateb.v15, stateb.q2[1], stateb.ϕ15)
# end

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
    ϕb = vrotate(ϕa + Δϕ_a, inv(qb) * qa)
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
# SET: Coordinates and Velocities
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
    # set_minimal_velocities!(pnode, cnode, joint; Δx=Δx, Δθ=Δθ, Δv=Δv, Δϕ=Δϕ)
    return nothing
end
