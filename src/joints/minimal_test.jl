
function child_velocities_alt(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    vb = (xb - xb1) / timestep
    ϕb = angular_velocity(qb1, qb, timestep)

    return [vb; ϕb]
end


function child_velocities_jacobian_velocity(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # Jacobians

    ∂vb∂Δv = ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra
    ∂vb∂Δϕ = ∂vrotate∂q(pb, qb1) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot
    ∂ϕb∂Δv = szeros(3, control_dimension(joint.translational))
    ∂ϕb∂Δϕ = ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot * timestep

    [
        ∂vb∂Δv ∂vb∂Δϕ;
        ∂ϕb∂Δv ∂ϕb∂Δϕ;
    ]
end

function ∂vb∂xa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)
   
    -1.0 / timestep * (I(3) + ∂vrotate∂p(pa + Atra * (minimal_coordinates(joint.translational, xa, qa, xb, qb) .- Δv * timestep), qa1) * Atra * minimal_coordinates_jacobian_configuration(:parent, joint.translational, xa, qa, xb, qb)[:, 1:3])
end

mech = get_slider()
initialize_slider!(mech, z1=0.1)
mech = get_pendulum()
initialize_pendulum!(mech, ϕ1=0.25 * π, ω1=0.1)
simulate!(mech, 0.5)
x = get_minimal_state(mech)

xa, va, qa, ϕa = current_configuration_velocity(mech.origin.state)
xb, qb = current_configuration(mech.bodies[1].state)

nu = control_dimension(mech.joints[1])
nu_tra = control_dimension(mech.joints[1].translational)
nu_rot = control_dimension(mech.joints[1].rotational)
child_velocities_alt(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    w, va, qa, ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), xa)

∂vb∂xa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    xa, w, qa, ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), va)

∂vb∂va(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    xa, va, UnitQuaternion(w..., false), ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), vector(qa))

∂vb∂qa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


function ∂vb∂va(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)
    1.0 * I(3)
end

function ∂vb∂qa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)

    J = -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:parent, joint.translational, xa, qa, xb, qb, attjac=false)[:, 3 .+ (1:4)]
    J += -1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
    J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2
    J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Lmat(qa1) * Rmat(qb * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Tmat()
    return J
end

∂vb∂qa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


function ∂vb∂ϕa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)

    # -1.0 / timestep * (vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1))
    # Lmat(qa) * Rmat(qoffset * Δq1) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2

    J = -1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2
    J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 
    
    return J
end

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    xa, va, qa, w, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), ϕa)[1:3, :]

∂vb∂ϕa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

function ∂ϕb∂xa(joint::JointConstraint{T},
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)

    return szeros(T, 3, 3)
end

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    w, va, qa, ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), xa)[3 .+ (1:3), :]

∂ϕb∂xa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


function ∂ϕb∂va(joint::JointConstraint{T},
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)

    return szeros(T, 3, 3)
end

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    xa, w, qa, ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), va)[3 .+ (1:3), :]

∂ϕb∂va(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

function ∂ϕb∂qa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)
    # ∂angular_velocity∂q1(qb1, qb, timestep) * qa1 * qoffset * Δq1

    J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2
    J += ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1) * Rmat(qb  * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Tmat() 
    return J
end

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    xa, va, UnitQuaternion(w..., false), ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), vector(qa))[3 .+ (1:3), :]

∂ϕb∂qa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


function ∂ϕb∂ϕa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
    Δq = inv(qoffset) * inv(qa) * qb 

    # step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)

    # step backward in time
    Δx1 = Δx .- Δv * timestep  
    Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
    qb1 = qa1 * qoffset * Δq1
    xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
    # finite-difference velocities
    # vb = (xb - xb1) / timestep
    # ϕb = angular_velocity(qb1, qb, timestep)

    J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 
    return J
end

FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
    xa, va, qa, w, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), ϕa)[3 .+ (1:3), :]

∂ϕb∂ϕa(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

child_velocities_jacobian_velocity(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# ∂vb∂Δv(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu.+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# ∂vb∂Δϕ(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu.+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# ∂ϕb∂Δv(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu.+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# ∂ϕb∂Δϕ(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu.+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# qa = rand(UnitQuaternion)
# qb = rand(UnitQuaternion)
# qc = rand(UnitQuaternion)
# p1 = rand(3)

# function test_rotation(x)
#     qd = qa * qb * qc * inv(axis_angle_to_quaternion(x))
#     vrotate(p1, qd)
# end

# r = rand(3)
# # r = vector(rand(UnitQuaternion))
# test_rotation(r)
# FiniteDiff.finite_difference_jacobian(test_rotation, r)
# ∂vrotate∂q(p1, qa * qb * qc * inv(axis_angle_to_quaternion(r))) * Lmat(qa * qb * qc) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(r)

