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
    Ta = impulse_transform_parent(joint, xa, qa, xb, qb)
    Tb = impulse_transform_child(joint, xa, qa, xb, qb)
    Faw = Ta[1:3,1:3] * Fτ
    τaa = Ta[4:6,1:3] * Fτ
    Fbw = Tb[1:3,1:3] * Fτ
    τbb = Tb[4:6,1:3] * Fτ
    return Faw, τaa, Fbw, τbb
end

@inline function input_jacobian_control_parent(joint::Translational, statea::State, stateb::State) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    # dFaw/dFτ
    # dτaa/dFτ
    Ta = impulse_transform_parent(joint, xa, qa, xb, qb)
    X = Ta[1:3,1:3]
    Q = 0.5*Ta[4:6,1:3]
    return [X; Q]
end

@inline function input_jacobian_control_child(joint::Translational, statea::State, stateb::State) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    # dFbw/dFτ
    # dτbb/dFτ
    Tb = impulse_transform_child(joint, xa, qa, xb, qb)
    X = Tb[1:3,1:3]
    Q = 0.5*Tb[4:6,1:3]
    return [X; Q]
end

@inline function input_jacobian_configuration_parent(joint::Translational{T}, statea::State, stateb::State) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
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
    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end

@inline function input_jacobian_configuration_child(joint::Translational{T}, statea::State, stateb::State) where T
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
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
    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
