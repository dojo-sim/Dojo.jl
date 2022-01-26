@inline function apply_input!(joint::Translational{T}, statea::State, stateb::State, timestep::T, clear::Bool) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)

    Faw, τaa, Fbw, τbb = apply_input(joint, joint.Fτ, xa, qa, xb, qb)
    statea.F2[end] += Faw
    statea.τ2[end] += τaa/2
    stateb.F2[end] += Fbw
    stateb.τ2[end] += τbb/2
    clear && (joint.Fτ = szeros(T,3))
    return
end

@inline function apply_input(joint::Translational{T}, F::AbstractVector, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    vertices = joint.vertices

    Faw = vrotate(-F, qa) # in the world frame
    Fbw = -Faw # in the world frame
    Faa = vrotate(Faw, inv(qa)) # in local frame
    Fbb = vrotate(Fbw, inv(qb)) # in local frame

    pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point in b frame
    ca_a = rotation_matrix(inv(qa)) * (xa) # body a com in a frame
    ra = pb_a - ca_a
    τaa = torque_from_force(Faa, ra) # in local coordinates
    τbb = torque_from_force(Fbb, vertices[2]) # in local coordinates
    return Faw, τaa, Fbw, τbb
end

@inline function input_jacobian_control_parent(joint::Translational, statea::State, stateb::State, timestep::T) where T
    vertices = joint.vertices
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)


    BFa = FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[1], joint.Fτ)
    Bτa = 0.5 * FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[2], joint.Fτ)

    return [BFa; Bτa]
end

@inline function input_jacobian_control_child(joint::Translational, statea::State, stateb::State, timestep::T) where T
    vertices = joint.vertices
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)

    BFb = FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[3], joint.Fτ)
    Bτb = 0.5 * FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[4], joint.Fτ)

    return [BFb; Bτb]
end

@inline function input_jacobian_configuration_parent(joint::Translational{T}, statea::State, stateb::State, timestep::T) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXa = FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[1], xa)
    FaQa = FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[1], [qa.w, qa.x, qa.y, qa.z])
    τaXa = 0.5 * FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[2], xa)
    τaQa = 0.5 * FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[2], [qa.w, qa.x, qa.y, qa.z])
    FbXa = FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[3], xa)
    FbQa = FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[3], [qa.w, qa.x, qa.y, qa.z])
    τbXa = 0.5 * FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[4], xa)
    τbQa = 0.5 * FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[4], [qa.w, qa.x, qa.y, qa.z])

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end

@inline function input_jacobian_configuration_child(joint::Translational{T}, statea::State, stateb::State, timestep::T) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[1], xb)
    FaQb = FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[1], [qb.w, qb.x, qb.y, qb.z])
    τaXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[2], xb)
    τaQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[2], [qb.w, qb.x, qb.y, qb.z])
    FbXb = FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[3], xb)
    FbQb = FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[3], [qb.w, qb.x, qb.y, qb.z])
    τbXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[4], xb)
    τbQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[4], [qb.w, qb.x, qb.y, qb.z])

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end