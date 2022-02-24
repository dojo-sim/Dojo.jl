################################################################################
# Coordinates
################################################################################
function minimal_coordinates(joint::JointConstraint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
	Δx = minimal_coordinates(joint.translational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
    Δθ = minimal_coordinates(joint.rotational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
	return [Δx; Δθ]
end

function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::JointConstraint, timestep;
        Δx::AbstractVector=szeros(control_dimension(joint.translational)),
        Δθ::AbstractVector=szeros(control_dimension(joint.rotational)))
    set_minimal_coordinates!(pnode, cnode, joint.rotational, timestep; Δθ=Δθ)
    set_minimal_coordinates!(pnode, cnode, joint.translational, timestep; Δx=Δx)
    return nothing
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

################################################################################
# Velocities
################################################################################
function minimal_velocities(joint::JointConstraint, pnode::Node, cnode::Node, timestep)
    Δv = minimal_velocities(joint.translational, pnode, cnode, timestep)
    Δϕ = minimal_velocities(joint.rotational, pnode, cnode, timestep)
    return [Δv; Δϕ]
end

function set_minimal_velocities!(pnode::Node, cnode::Node, joint::JointConstraint, timestep;
        Δv=szeros(control_dimension(joint.translational)),
        Δϕ=szeros(control_dimension(joint.rotational)))
	vb, ϕb = set_child_velocities(joint, 
        initial_configuration_velocity(pnode.state)...,
	 	current_configuration(cnode.state)..., timestep, Δv=Δv, Δϕ=Δϕ)
	set_velocity!(cnode; v=vb, ω=ϕb)
	set_previous_configuration!(cnode, timestep)
    return nothing
end

function set_child_velocities(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arotᵀ = zerodimstaticadjoint(nullspace_mask(rot))
    Atraᵀ = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb

    # 1 step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # Finite difference
    Δx1 = Δx .- Δv * timestep
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arotᵀ * Δϕ * timestep))

    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atraᵀ*Δx1, qa1) - vrotate(pb, qb1)
    
    # Finite difference
    vb = (xb - xb1) / timestep
    ϕb = angular_velocity(qb1, qb, timestep)
    return vb, ϕb
end

function set_velocity!(mechanism, joint::JointConstraint{T,N,Nc}, vϕ) where {T,N,Nc}
    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    # unpack
    Δv = vϕ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δϕ = vϕ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    
    set_minimal_velocities!(body1, body2, joint, mechanism.timestep, Δv=Δv, Δϕ=Δϕ)
    return body2.state.v15, body2.state.ϕ15
end

################################################################################
# Coordinates and Velocities
################################################################################
function minimal_coordinates_velocities(joint::JointConstraint,
    pnode::Node, cnode::Node, timestep)
    Δxθ = minimal_coordinates(joint, pnode, cnode)
    Δvϕ = minimal_velocities(joint, pnode, cnode, timestep)
    return [Δxθ; Δvϕ]
end

function set_minimal_coordinates_velocities!(mechanism::Mechanism, joint::JointConstraint;
        xmin::AbstractVector=szeros(2*control_dimension(joint)))
    pnode = get_body(mechanism, joint.parent_id)
    cnode = get_body(mechanism, joint.child_id)
    xa = pnode.state.x2[1] 
    va = pnode.state.v15 
    qa = pnode.state.q2[1] 
    ϕa = pnode.state.ϕ15 
    za = [xa; va; vector(qa); ϕa]

    set_minimal_coordinates_velocities!(pnode, cnode, joint, mechanism.timestep; xmin=xmin, zp=za)
end

function set_minimal_coordinates_velocities!(pbody::Node, cbody::Node, joint::JointConstraint, timestep;
        zp::AbstractVector=[pbody.state.x2[1]; pbody.state.v15; vector(pbody.state.q2[1]); pbody.state.ϕ15],
        xmin::AbstractVector=szeros(2*control_dimension(joint)))
    
    nu = control_dimension(joint)

    Δx = xmin[SUnitRange(joint.minimal_index[1]...)]
    Δθ = xmin[SUnitRange(joint.minimal_index[2]...)]

    Δv = xmin[nu .+ SUnitRange(joint.minimal_index[1]...)]
    Δϕ = xmin[nu .+ SUnitRange(joint.minimal_index[2]...)]

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    # parent state
    xa = SVector{3}(zp[1:3])
    va = SVector{3}(zp[3 .+ (1:3)])
    qa = UnitQuaternion(zp[6 .+ (1:4)]..., false)
    ϕa = SVector{3}(zp[10 .+ (1:3)])

    # positions
    Δq = axis_angle_to_quaternion(Arot * Δθ)
	qb = qa * qoffset * Δq
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    # previous configuration
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # finite-difference configuration
    Δx1 = Δx .- Δv * timestep
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))

    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocity
    vb = (xb - xb1) / timestep
    ϕb = angular_velocity(qb1, qb, timestep)
    
    # set child state
    cbody.state.x2[1] = xb 
    cbody.state.v15 = vb 
    cbody.state.q2[1] = qb 
    cbody.state.ϕ15 = ϕb 

    return [xb; vb; vector(qb); ϕb]
end

function minimal_coordinates_velocities_jacobian_parent(joint::JointConstraint{T},
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    # positions
    Δq = axis_angle_to_quaternion(Arot * Δθ)
	qb = qa * qoffset * Δq
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # Jacobians

    ∂xb∂xa = 1.0 * I(3)

    ∂xb∂va = szeros(T, 3, 3)

    ∂xb∂qa = ∂vrotate∂q(pa + Atra * Δx, qa)
    ∂xb∂qa += -∂vrotate∂q(pb, qb) * Rmat(qoffset * axis_angle_to_quaternion(Arot * Δθ))

    ∂xb∂ϕa = szeros(T, 3, 3)

    ∂qb∂xa = szeros(T, 4, 3)

    ∂qb∂va = szeros(T, 4, 3)

    ∂qb∂qa = Rmat(qoffset * axis_angle_to_quaternion(Arot * Δθ))

    ∂qb∂ϕa = szeros(T, 4, 3)

    ∂vb∂xa = szeros(T, 3, 3)
    
    ∂vb∂va = 1.0 * I(3)

    ∂vb∂qa = 1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx, qa)
    ∂vb∂qa -= 1.0 / timestep * ∂vrotate∂q(pb, qb) * Rmat(qoffset * Δq)
    ∂vb∂qa += -1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
    ∂vb∂qa += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2 

    ∂vb∂ϕa = 1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2
    ∂vb∂ϕa += -1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep)  * timestep / 2
    
    ∂ϕb∂xa = szeros(T, 3, 3) 

    ∂ϕb∂va = szeros(T, 3, 3)

    ∂ϕb∂qa = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2
    ∂ϕb∂qa += ∂angular_velocity∂q2(qb1, qb, timestep) * Rmat(qoffset * Δq)
    
    ∂ϕb∂ϕa = -1.0 * ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 

    [
        ∂xb∂xa ∂xb∂va ∂xb∂qa ∂xb∂ϕa;
        ∂vb∂xa ∂vb∂va ∂vb∂qa ∂vb∂ϕa;
        ∂qb∂xa ∂qb∂va ∂qb∂qa ∂qb∂ϕa;
        ∂ϕb∂xa ∂ϕb∂va ∂ϕb∂qa ∂ϕb∂ϕa;
    ]
end

function minimal_coordinates_velocities_jacobian_minimal(joint::JointConstraint{T},
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    # positions
    Δq = axis_angle_to_quaternion(Arot * Δθ)
	qb = qa * qoffset * Δq
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)
    
    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δq = axis_angle_to_quaternion(Arot * Δθ)
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # Jacobians

    ∂xb∂Δx = ∂vrotate∂p(pa + Atra * Δx, qa) * Atra
    ∂xb∂Δθ = -∂vrotate∂q(pb, qb) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    ∂xb∂Δv = szeros(T, 3, control_dimension(joint.translational))
    ∂xb∂Δϕ = szeros(T, 3, control_dimension(joint.rotational))

    ∂qb∂Δx = szeros(T, 4, control_dimension(joint.translational))
    ∂qb∂Δθ = Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    ∂qb∂Δv = szeros(T, 4, control_dimension(joint.translational)) 
    ∂qb∂Δϕ = szeros(T, 4, control_dimension(joint.rotational))

    ∂vb∂Δx = 1 / timestep * ∂vrotate∂p(pa + Atra * Δx, qa) * Atra
    ∂vb∂Δx += -1 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra
    
    ∂vb∂Δθ = -1.0 / timestep * ∂vrotate∂q(pb, qb) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    ∂vb∂Δθ += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    
    ∂vb∂Δv = ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra
    
    ∂vb∂Δϕ = ∂vrotate∂q(pb, qb1) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot
    
    ∂ϕb∂Δx = szeros(T, 3, control_dimension(joint.translational))

    ∂ϕb∂Δθ = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    ∂ϕb∂Δθ += ∂angular_velocity∂q2(qb1, qb, timestep) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    
    ∂ϕb∂Δv = szeros(3, control_dimension(joint.translational))
    
    ∂ϕb∂Δϕ = ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot * timestep

    [
        ∂xb∂Δx ∂xb∂Δθ ∂xb∂Δv ∂xb∂Δϕ;
        ∂vb∂Δx ∂vb∂Δθ ∂vb∂Δv ∂vb∂Δϕ;
        ∂qb∂Δx ∂qb∂Δθ ∂qb∂Δv ∂qb∂Δϕ;
        ∂ϕb∂Δx ∂ϕb∂Δθ ∂ϕb∂Δv ∂ϕb∂Δϕ;
    ]
end

function set_minimal_coordinates_velocities_jacobian_minimal(pnode::Node, cnode::Node,
    joint::JointConstraint, timestep)

    xmin = minimal_coordinates_velocities(joint, pnode, cnode, timestep)
    qb = current_orientation(cnode.state)
    xa, va, qa, ϕa = initial_configuration_velocity(pnode.state)
    za = [xa; va; vector(qa); ϕa]

    nu = control_dimension(joint)
    Δx = xmin[SUnitRange(joint.minimal_index[1]...)]
    Δθ = xmin[SUnitRange(joint.minimal_index[2]...)]
    Δv = xmin[nu .+ SUnitRange(joint.minimal_index[1]...)]
    Δϕ = xmin[nu .+ SUnitRange(joint.minimal_index[2]...)]
    J = minimal_coordinates_velocities_jacobian_minimal(joint,
                xa, va, qa, ϕa,
                timestep;
                Δx=Δx,
                Δθ=Δθ,
                Δv=Δv,
                Δϕ=Δϕ)
    J = cat(Diagonal(sones(6)), LVᵀmat(qb)', Diagonal(sones(3)), dims=(1,2)) * J
    return J
end

function set_minimal_coordinates_velocities_jacobian_parent(pnode::Node, cnode::Node,
    joint::JointConstraint, timestep)
    xmin = minimal_coordinates_velocities(joint, pnode, cnode, timestep)
    qb = current_orientation(cnode.state)
    xa, va, qa, ϕa = initial_configuration_velocity(pnode.state)
    za = [xa; va; vector(qa); ϕa]

    nu = control_dimension(joint)
    Δx = xmin[SUnitRange(joint.minimal_index[1]...)]
    Δθ = xmin[SUnitRange(joint.minimal_index[2]...)]
    Δv = xmin[nu .+ SUnitRange(joint.minimal_index[1]...)]
    Δϕ = xmin[nu .+ SUnitRange(joint.minimal_index[2]...)]
    J = minimal_coordinates_velocities_jacobian_parent(joint,
                xa, va, qa, ϕa,
                timestep;
                Δx=Δx,
                Δθ=Δθ,
                Δv=Δv,
                Δϕ=Δϕ)
    J = cat(Diagonal(sones(6)), LVᵀmat(qb)', Diagonal(sones(3)), dims=(1,2)) * J
    J = J * cat(Diagonal(sones(6)), LVᵀmat(qa), Diagonal(sones(3)), dims=(1,2))
    return J
end