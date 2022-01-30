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

@inline function apply_input(joint::Translational{T}, Fτ::AbstractVector,
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    # vertices = joint.vertices
    # Faw = vrotate(-Fτ, qa) # in the world frame
    # Fbw = -Faw # in the world frame
    # Faa = vrotate(Faw, inv(qa)) # in local frame
    # Fbb = vrotate(Fbw, inv(qb)) # in local frame
    #
    # pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point in b frame
    # ca_a = rotation_matrix(inv(qa)) * (xa) # body a com in a frame
    # ra = pb_a - ca_a
    # τaa = torque_from_force(Faa, ra) # in local coordinates
    # τbb = torque_from_force(Fbb, vertices[2]) # in local coordinates

    Ta = impulse_transform_parent(joint, xa, qa, xb, qb)
    Tb = impulse_transform_child(joint, xa, qa, xb, qb)
    Faw = Ta[1:3,1:3] * Fτ
    τaa = Ta[4:6,1:3] * Fτ
    Fbw = Tb[1:3,1:3] * Fτ
    τbb = Tb[4:6,1:3] * Fτ
    # @show norm(Faw-Faw0)
    # @show norm(τaa-τaa0)
    # @show norm(Fbw-Fbw0)
    # @show norm(τbb-τbb0)
    return Faw, τaa, Fbw, τbb
end

@inline function input_jacobian_control_parent(joint::Translational, statea::State, stateb::State, timestep::T) where T
    # vertices = joint.vertices
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    # BFa = FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[1], joint.Fτ)
    # Bτa = 0.5 * FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[2], joint.Fτ)

    # dFaw/dFτ
    # dτaa/dFτ
    Ta = impulse_transform_parent(joint, xa, qa, xb, qb)
    X = Ta[1:3,1:3]
    Q = 0.5*Ta[4:6,1:3]
    # @show norm(BFa - BFa0)
    # @show norm(Bτa - Bτa0)
    return [X; Q]
end

@inline function input_jacobian_control_child(joint::Translational, statea::State, stateb::State, timestep::T) where T
    # vertices = joint.vertices
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    #
    # BFb = FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[3], joint.Fτ)
    # Bτb = 0.5 * FiniteDiff.finite_difference_jacobian(F -> apply_input(joint, F, xa, qa, xb, qb)[4], joint.Fτ)

    # dFbw/dFτ
    # dτbb/dFτ
    Tb = impulse_transform_child(joint, xa, qa, xb, qb)
    X = Tb[1:3,1:3]
    Q = 0.5*Tb[4:6,1:3]
    # @show norm(BFb - BFb0)
    # @show norm(Bτb - Bτb0)
    return [X; Q]
end

@inline function input_jacobian_configuration_parent(joint::Translational{T}, statea::State, stateb::State, timestep::T) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    # F = joint.Fτ
    # vertices = joint.vertices

    # FaXa = FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[1], xa)
    # FaQa = FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[1], [qa.w, qa.x, qa.y, qa.z])
    # τaXa = 0.5 * FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[2], xa)
    # τaQa = 0.5 * FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[2], [qa.w, qa.x, qa.y, qa.z])
    # FbXa = FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[3], xa)
    # FbQa = FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[3], [qa.w, qa.x, qa.y, qa.z])
    # τbXa = 0.5 * FiniteDiff.finite_difference_jacobian(xa -> apply_input(joint, F, xa, qa, xb, qb)[4], xa)
    # τbQa = 0.5 * FiniteDiff.finite_difference_jacobian(qa -> apply_input(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[4], [qa.w, qa.x, qa.y, qa.z])

    # d[Faw;2τaa]/d[xa,qa]
    ∇aa = impulse_transform_parent_jacobian_parent(joint, xa, qa, xb, qb, joint.Fτ)
    FaXa = ∇aa[1:3,1:3]
    FaQa = ∇aa[1:3,4:6]
    τaXa = 0.5*∇aa[4:6,1:3]
    τaQa = 0.5*∇aa[4:6,4:6]
    # d[Fbw;2τbb]/d[xa,qa]
    ∇ba = impulse_transform_child_jacobian_parent(joint, xa, qa, xb, qb, joint.Fτ)
    FbXa = ∇ba[1:3,1:3]
    FbQa = ∇ba[1:3,4:6]
    τbXa = 0.5*∇ba[4:6,1:3]
    τbQa = 0.5*∇ba[4:6,4:6]

    # @show norm(FaXa0 - FaXa, Inf)
    # @show norm(FaQa0 - FaQa*LVᵀmat(qa), Inf)
    # @show norm(τaXa0 - τaXa, Inf)
    # @show norm(τaQa0 - τaQa*LVᵀmat(qa), Inf)
    # @show norm(FbXa0 - FbXa, Inf)
    # @show norm(FbQa0 - FbQa*LVᵀmat(qa), Inf)
    # @show norm(τbXa0 - τbXa, Inf)
    # @show norm(τbQa0 - τbQa*LVᵀmat(qa), Inf)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end

@inline function input_jacobian_configuration_child(joint::Translational{T}, statea::State, stateb::State, timestep::T) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    # F = joint.Fτ
    # vertices = joint.vertices

    # FaXb = FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[1], xb)
    # FaQb = FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[1], [qb.w, qb.x, qb.y, qb.z])
    # τaXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[2], xb)
    # τaQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[2], [qb.w, qb.x, qb.y, qb.z])
    # FbXb = FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[3], xb)
    # FbQb = FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[3], [qb.w, qb.x, qb.y, qb.z])
    # τbXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> apply_input(joint, F, xa, qa, xb, qb)[4], xb)
    # τbQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> apply_input(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[4], [qb.w, qb.x, qb.y, qb.z])


    # d[Faw;2τaa]/d[xb,qb]
    ∇ab = impulse_transform_parent_jacobian_child(joint, xa, qa, xb, qb, joint.Fτ)
    FaXb = ∇ab[1:3,1:3]
    FaQb = ∇ab[1:3,4:6]
    τaXb = 0.5*∇ab[4:6,1:3]
    τaQb = 0.5*∇ab[4:6,4:6]
    # d[Fbw;2τbb]/d[xb,qb]
    ∇bb = impulse_transform_child_jacobian_child(joint, xa, qa, xb, qb, joint.Fτ)
    FbXb = ∇bb[1:3,1:3]
    FbQb = ∇bb[1:3,4:6]
    τbXb = 0.5*∇bb[4:6,1:3]
    τbQb = 0.5*∇bb[4:6,4:6]

    # @show norm(FaXb0 - FaXb, Inf)
    # @show norm(FaQb0 - FaQb*LVᵀmat(qb), Inf)
    # @show norm(τaXb0 - τaXb, Inf)
    # @show norm(τaQb0 - τaQb*LVᵀmat(qb), Inf)
    # @show norm(FbXb0 - FbXb, Inf)
    # @show norm(FbQb0 - FbQb*LVᵀmat(qb), Inf)
    # @show norm(τbXb0 - τbXb, Inf)
    # @show norm(τbQb0 - τbQb*LVᵀmat(qb), Inf)
    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
