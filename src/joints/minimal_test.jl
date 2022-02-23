
# function child_velocities_alt(joint::JointConstraint,
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational)))

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     vb = (xb - xb1) / timestep
#     ϕb = angular_velocity(qb1, qb, timestep)

#     return [vb; ϕb]
# end

function child_velocities_jacobian_child_configurations(joint::JointConstraint{T},
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
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
    ∂vb∂xb = 1.0 / timestep * I(3)
    ∂vb∂xb += -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:child, joint.translational, xa, qa, xb, qb)[:, 1:3]
    
    ∂vb∂qb = -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:child, joint.translational, xa, qa, xb, qb, attjac=false)[:, 3 .+ (1:4)]
    ∂vb∂qb += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset * inv(qoffset) * inv(qa))

    ∂ϕb∂xb = szeros(T, 3, 3)

    ∂ϕb∂qb = ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1 * qoffset * inv(qoffset) * inv(qa)) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep)))
    ∂ϕb∂qb += ∂angular_velocity∂q2(qb1, qb, timestep)

    [
        ∂vb∂xb ∂vb∂qb;
        ∂ϕb∂xb ∂ϕb∂qb;
    ]
end

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    w[1:3], UnitQuaternion(w[3 .+ (1:4)]..., false),
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))...), [xb; vector(qb)])[:, :]

child_velocities_jacobian_child_configurations(mech.joints[1], 
    xa, va, qa, ϕa,
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


function ∂vb∂xb(joint::JointConstraint,
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

    # return [vb; ϕb]
    J = 1.0 / timestep * I(3)
    J += -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:child, joint.translational, xa, qa, xb, qb)[:, 1:3]
    return J
end

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    w, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[1:nu_tra]),
    Δϕ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), xb)[1:3, :]

∂vb∂xb(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

function ∂vb∂qb(joint::JointConstraint,
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

    # return [vb; ϕb]
    -1.0 / timestep * (vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1))

    J = -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:child, joint.translational, xa, qa, xb, qb, attjac=false)[:, 3 .+ (1:4)]
    J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset * inv(qoffset) * inv(qa))
    return J
end

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    xb, UnitQuaternion(w..., false),
    mech.timestep,
    Δv=SVector{nu_tra}(x[1:nu_tra]),
    Δϕ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), vector(qb))[1:3, :]
    ∂∂
∂vb∂qb(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

function ∂ϕb∂xb(joint::JointConstraint{T},
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

    # return [vb; ϕb]

    return szeros(T, 3, 3)
end

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    w, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[1:nu_tra]),
    Δϕ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), xb)[4:6, :]

∂ϕb∂xb(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))
   

function ∂ϕb∂qb(joint::JointConstraint{T},
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

    # return [vb; ϕb]
    J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset * inv(qoffset) * inv(qa))
    J += ∂angular_velocity∂q2(qb1, qb, timestep)
    return J
end

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    xb, UnitQuaternion(w..., false),
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu .+ (nu_tra .+ (1:nu_rot))]))...), vector(qb))[4:6, :]

∂ϕb∂qb(mech.joints[1], 
    xa, va, qa, ϕa,
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))
   
FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    w[1:3], UnitQuaternion(w[3 .+ (1:4)]..., false),
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu .+ (nu_tra .+ (1:nu_rot))]))...), [xb; vector(qb)])[4:6, :]

∂ϕb∂qb(mech.joints[1], 
    xa, va, qa, ϕa,
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))
       
    

function ∂vb∂qb(joint::JointConstraint,
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

    # return [vb; ϕb]
    1.0 / timestep * I(3)
end

function ∂vb∂qa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    J = 1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx, qa)
    J -= 1.0 / timestep * ∂vrotate∂q(pb, qb) * Rmat(qoffset * Δq)
    J += -1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
    J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2 

    return J
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=xmin, zp=w), za)[collect(3 .+ (1:3)), 6 .+ (1:4)]

∂vb∂qa(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)


function ∂vb∂ϕa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    -1.0 / timestep * (vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1))

    J = 1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2
    J += -1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep)  * timestep / 2
    return J
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=xmin, zp=w), za)[collect(3 .+ (1:3)), 10 .+ (1:3)]

∂vb∂ϕa(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)


function ∂ϕb∂qa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2
    J += ∂angular_velocity∂q2(qb1, qb, timestep) * Rmat(qoffset * Δq)
    return J
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=xmin, zp=w), za)[collect(10 .+ (1:3)), 6 .+ (1:4)]

∂ϕb∂qa(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)

function ∂ϕb∂ϕa(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    J = -1.0 * ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2
    return J
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=xmin, zp=w), za)[collect(10 .+ (1:3)), 10 .+ (1:3)]

∂ϕb∂ϕa(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)

function ∂vb∂Δx(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    1 / timestep * (vrotate(pa + Atra * Δx, qa) - vrotate(pa + Atra * Δx1, qa1))

    J = 1 / timestep * ∂vrotate∂p(pa + Atra * Δx, qa) * Atra
    J += -1 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra

    return J 
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=w, zp=za), xmin)[collect(3 .+ (1:3)), 1:nu_tra]

∂vb∂Δx(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)

function ∂vb∂Δθ(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    J = -1.0 / timestep * ∂vrotate∂q(pb, qb) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    
    1 / timestep * (- vrotate(pb, qb) - (- vrotate(pb, qb1)))
    return J 
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=w, zp=za), xmin)[collect(3 .+ (1:3)), nu_tra .+ (1:nu_rot)]

∂vb∂Δθ(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)

function ∂ϕb∂Δθ(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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

    J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    J += ∂angular_velocity∂q2(qb1, qb, timestep) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    return J 
end

# velocity wrt maximal
FiniteDiff.finite_difference_jacobian(w -> set_minimal_coordinates_velocities!(bodya, bodyb, joint, timestep, xmin=w, zp=za), xmin)[collect(10 .+ (1:3)), nu_tra .+ (1:nu_rot)]

∂ϕb∂Δθ(joint,
    xa, va, qa, ϕa,
    xb, qb, 
    timestep;
    Δx=Δx,
    Δθ=Δθ,
    Δv=Δv,
    Δϕ=Δϕ)
    


function ∂∂(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

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


    ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot * timestep
end

# mech = get_slider()
# initialize_slider!(mech, z1=0.1)
mech = get_pendulum()
initialize_pendulum!(mech, ϕ1=0.25 * π, ω1=0.1)
# mech = get_halfcheetah()
# initialize_halfcheetah!(mech)
simulate!(mech, 0.5)
x = get_minimal_state(mech)

xa, va, qa, ϕa = current_configuration_velocity(mech.origin.state)
xb, qb = current_configuration(mech.bodies[1].state)

nu = control_dimension(mech.joints[1])
nu_tra = control_dimension(mech.joints[1].translational)
nu_rot = control_dimension(mech.joints[1].rotational)
vcat(set_child_velocities(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))...)

set_child_velocities_jacobian_velocity(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    xa, va, qa, ϕa, 
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(w[1:nu_tra]),
    Δϕ=SVector{nu_rot}(w[nu_tra .+ (1:nu_rot)]))...), x[nu .+ (1:nu)])

set_child_velocities_jacobian_parent(mech.joints[1], 
    current_configuration_velocity(mech.origin.state)...,
    current_configuration(mech.bodies[1].state)...,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_velocities(mech.joints[1], 
    SVector{3}(w[1:3]), SVector{3}(w[3 .+ (1:3)]), UnitQuaternion(w[6 .+ (1:4)]..., false), SVector{3}(w[10 .+ (1:3)]),
    xb, qb,
    mech.timestep,
    Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
    Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))...), [xa; va; vector(qa); ϕa])


a = 1


# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     w, va, qa, ϕa, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), xa)

# ∂vb∂xa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     xa, w, qa, ϕa, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), va)

# ∂vb∂va(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     xa, va, UnitQuaternion(w..., false), ϕa, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), vector(qa))

# ∂vb∂qa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


# function ∂vb∂va(joint::JointConstraint,
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational)))

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)
#     1.0 * I(3)
# end

# function ∂vb∂qa(joint::JointConstraint,
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational)))

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)

#     J = -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:parent, joint.translational, xa, qa, xb, qb, attjac=false)[:, 3 .+ (1:4)]
#     J += -1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
#     J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2
#     J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Lmat(qa1) * Rmat(qb * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Tmat()
#     return J
# end

# ∂vb∂qa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


# function ∂vb∂ϕa(joint::JointConstraint,
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational)))

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)

#     # -1.0 / timestep * (vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1))
#     # Lmat(qa) * Rmat(qoffset * Δq1) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2

#     J = -1.0 / timestep * ∂vrotate∂q(pa + Atra * Δx1, qa1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2
#     J += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 
    
#     return J
# end

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     xa, va, qa, w, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), ϕa)[1:3, :]

# ∂vb∂ϕa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# function ∂ϕb∂xa(joint::JointConstraint{T},
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational))) where T

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)

#     return szeros(T, 3, 3)
# end

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     w, va, qa, ϕa, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), xa)[3 .+ (1:3), :]

# ∂ϕb∂xa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


# function ∂ϕb∂va(joint::JointConstraint{T},
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational))) where T

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)

#     return szeros(T, 3, 3)
# end

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     xa, w, qa, ϕa, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), va)[3 .+ (1:3), :]

# ∂ϕb∂va(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))

# function ∂ϕb∂qa(joint::JointConstraint,
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational)))

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)
#     # ∂angular_velocity∂q1(qb1, qb, timestep) * qa1 * qoffset * Δq1

#     J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(quaternion_map(-ϕa, timestep) * qoffset * Δq1) * timestep / 2
#     J += ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1) * Rmat(qb  * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Tmat() 
#     return J
# end

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     xa, va, UnitQuaternion(w..., false), ϕa, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), vector(qa))[3 .+ (1:3), :]

# ∂ϕb∂qa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


# function ∂ϕb∂ϕa(joint::JointConstraint,
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, timestep;
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational)))

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)
#     Δq = inv(qoffset) * inv(qa) * qb 

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     # vb = (xb - xb1) / timestep
#     # ϕb = angular_velocity(qb1, qb, timestep)

#     J = ∂angular_velocity∂q1(qb1, qb, timestep) * Rmat(qoffset * Δq1) * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 
#     return J
# end

# FiniteDiff.finite_difference_jacobian(w -> child_velocities_alt(mech.joints[1], 
#     xa, va, qa, w, 
#     xb, qb,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)])), ϕa)[3 .+ (1:3), :]

# ∂ϕb∂ϕa(mech.joints[1], 
#     current_configuration_velocity(mech.origin.state)...,
#     current_configuration(mech.bodies[1].state)...,
#     mech.timestep,
#     Δv=SVector{nu_tra}(x[nu .+ (1:nu_tra)]),
#     Δϕ=SVector{nu_rot}(x[nu + nu_tra .+ (1:nu_rot)]))


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




# FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
#     xa, qa, 
#     mech.timestep,
#     Δx=SVector{nu_tra}(x[1:nu_tra]),
#     Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...))[4:6, :]



set_child_configurations(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

child_configurations_jacobian_configuration(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    w[1:3], UnitQuaternion(w[3 .+ (1:4)]..., false), 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), [xa; vector(qa)])



child_configurations_jacobian_minimal(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(w[1:nu_tra]),
    Δθ=SVector{nu_rot}(w[nu_tra .+ (1:nu_rot)]))...), x[1:nu])

function ∂xb∂xa(joint::JointConstraint,
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return 1.0 * I(3)
end

∂xb∂xa(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    w, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), xa)[1:3, :]

function ∂xb∂qa(joint::JointConstraint,
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    J = ∂vrotate∂q(pa + Atra * Δx, qa)
    J += -∂vrotate∂q(pb, qb) * Rmat(qoffset * axis_angle_to_quaternion(Arot * Δθ))
    return J
end

∂xb∂qa(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, UnitQuaternion(w..., false), 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), vector(qa))[1:3, :]

function ∂qb∂xa(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return szeros(T, 4, 3)
end

∂qb∂xa(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    w, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), vector(xa))[3 .+ (1:4), :]

function ∂qb∂qa(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return Rmat(qoffset * axis_angle_to_quaternion(Arot * Δθ))
end
    
∂qb∂qa(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, UnitQuaternion(w..., false), 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), vector(qa))[3 .+ (1:4), :]
    

function ∂xb∂Δx(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return ∂vrotate∂p(pa + Atra * Δx, qa) * Atra
end
    
∂xb∂Δx(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=w,
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), x[1:nu_tra])[1:3, :]
       
function ∂xb∂Δθ(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return -∂vrotate∂q(pb, qb) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
end
    
∂xb∂Δθ(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(w))...), x[nu_tra .+ (1:nu_rot)])[1:3, :]
     
    
function ∂qb∂Δx(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return szeros(T, 4, 3)
end
    
∂qb∂Δx(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(w),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))...), x[1:nu_tra])[3 .+ (1:4), :]

function ∂qb∂Δθ(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
end
    
∂qb∂Δθ(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(x[nu_tra .+ (1:nu_rot)]))

FiniteDiff.finite_difference_jacobian(w -> vcat(set_child_configurations(mech.joints[1], 
    xa, qa, 
    mech.timestep,
    Δx=SVector{nu_tra}(x[1:nu_tra]),
    Δθ=SVector{nu_rot}(w))...), x[nu_tra .+ (1:nu_rot)])[3 .+ (1:4), :]
             
############
# function child_velocities_jacobian_child_configurations(joint::JointConstraint{T},
#     xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#     xb::AbstractVector, qb::UnitQuaternion, 
#     timestep;
#     Δx=szeros(control_dimension(joint.translational)),
#     Δθ=szeros(control_dimension(joint.rotational)),
#     Δv=szeros(control_dimension(joint.translational)),
#     Δϕ=szeros(control_dimension(joint.rotational))) where T

#     rot = joint.rotational
#     tra = joint.translational
#     pa = tra.vertices[1]
#     pb = tra.vertices[2]
#     qoffset = rot.qoffset
#     Arot = zerodimstaticadjoint(nullspace_mask(rot))
#     Atra = zerodimstaticadjoint(nullspace_mask(tra))

#     # step backward in time
#     xa1 = next_position(xa, -va, timestep)
#     qa1 = next_orientation(qa, -ϕa, timestep)

#     # step backward in time
#     Δq = axis_angle_to_quaternion(Arot * Δθ)
#     Δx1 = Δx .- Δv * timestep  
#     Δq1 = Δq * inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))
    
#     qb1 = qa1 * qoffset * Δq1
#     xb1 = xa1 + vrotate(pa + Atra * Δx1, qa1) - vrotate(pb, qb1)
    
#     # finite-difference velocities
#     ∂vb∂xb = 1.0 / timestep * I(3)
#     ∂vb∂xb += -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:child, joint.translational, xa, qa, xb, qb)[:, 1:3]
    
#     ∂vb∂qb = -1.0 / timestep * ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra * minimal_coordinates_jacobian_configuration(:child, joint.translational, xa, qa, xb, qb, attjac=false)[:, 3 .+ (1:4)]
#     ∂vb∂qb += 1.0 / timestep * ∂vrotate∂q(pb, qb1) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep))) * Lmat(qa1 * qoffset * inv(qoffset) * inv(qa))

#     ∂ϕb∂xb = szeros(T, 3, 3)

#     ∂ϕb∂qb = ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1 * qoffset * inv(qoffset) * inv(qa)) * Rmat(inv(axis_angle_to_quaternion(Arot * Δϕ * timestep)))
#     ∂ϕb∂qb += ∂angular_velocity∂q2(qb1, qb, timestep)

#     [
#         ∂vb∂xb ∂vb∂qb;
#         ∂ϕb∂xb ∂ϕb∂qb;
#     ]
# end

function child_velocities_jacobian_minimal(joint::JointConstraint,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)),
    Δv=szeros(control_dimension(joint.translational)),
    Δϕ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

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

    ∂vb∂Δv = ∂vrotate∂p(pa + Atra * Δx1, qa1) * Atra
    ∂vb∂Δϕ = ∂vrotate∂q(pb, qb1) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot
    ∂ϕb∂Δv = szeros(3, control_dimension(joint.translational))
    ∂ϕb∂Δϕ = ∂angular_velocity∂q1(qb1, qb, timestep) * Lmat(qa1 * qoffset * Δq) * Tmat() * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δϕ * timestep) * Arot * timestep

    [
        ∂vb∂Δv ∂vb∂Δϕ;
        ∂ϕb∂Δv ∂ϕb∂Δϕ;
    ]
end

function child_velocities_jacobian_parent_configuration(joint::JointConstraint{T},
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, qb::UnitQuaternion, timestep;
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
        ∂vb∂xa ∂vb∂va ∂vb∂qa ∂vb∂ϕa;
        ∂ϕb∂xa ∂ϕb∂va ∂ϕb∂qa ∂ϕb∂ϕa;
    ]
end

function set_child_configurations(joint::JointConstraint,
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational)))

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

	qb = qa * qoffset * axis_angle_to_quaternion(Arot * Δθ)
    xb = xa + vrotate(pa + Atra * Δx, qa) - vrotate(pb, qb)

    return xb, vector(qb)
end

function child_configurations_jacobian_parent_configuration(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    ∂xb∂xa = 1.0 * I(3)

    ∂xb∂qa = ∂vrotate∂q(pa + Atra * Δx, qa)
    ∂xb∂qa += -∂vrotate∂q(pb, qb) * Rmat(qoffset * axis_angle_to_quaternion(Arot * Δθ))

    ∂qb∂xa = szeros(T, 4, 3)
    ∂qb∂qa = Rmat(qoffset * axis_angle_to_quaternion(Arot * Δθ))

	[
        ∂xb∂xa ∂xb∂qa;
        ∂qb∂xa ∂qb∂qa;
    ]
end

function child_configurations_jacobian_minimal(joint::JointConstraint{T},
    xa::AbstractVector, qa::UnitQuaternion,
    timestep;
    Δx=szeros(control_dimension(joint.translational)),
    Δθ=szeros(control_dimension(joint.rotational))) where T

    rot = joint.rotational
    tra = joint.translational
    pa = tra.vertices[1]
    pb = tra.vertices[2]
    qoffset = rot.qoffset
    Arot = zerodimstaticadjoint(nullspace_mask(rot))
    Atra = zerodimstaticadjoint(nullspace_mask(tra))

    ∂xb∂Δx = ∂vrotate∂p(pa + Atra * Δx, qa) * Atra
    ∂xb∂Δθ = -∂vrotate∂q(pb, qb) * Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot
    ∂qb∂Δx = szeros(T, 4, control_dimension(joint.translational))
    ∂qb∂Δθ = Lmat(qa * qoffset) * ∂axis_angle_to_quaternion∂axis_angle(Arot * Δθ) * Arot

    [
        ∂xb∂Δx ∂xb∂Δθ;
        ∂qb∂Δx ∂qb∂Δθ;
    ]
end
