@inline function apply_input!(joint::Rotational{T}, statea::State, stateb::State, timestep::T, clear::Bool) where T
    τ = joint.Fτ
    _, qa = current_configuration(statea)
    _, qb = current_configuration(stateb)

    statea.τ2[end] += -τ
    stateb.τ2[end] += vrotate(vrotate(τ, qa),inv(qb))
    clear && (joint.Fτ = szeros(T,3))
    return
end

@inline function input_jacobian_control_parent(joint::Rotational{T}, statea::State, stateb::State) where T
    BFa = szeros(T, 3, 3)
    Bτa = -I
    return [BFa; Bτa]
end

@inline function input_jacobian_control_child(joint::Rotational{T}, statea::State, stateb::State) where T
    _, qa = current_configuration(statea)
    _, qb = current_configuration(stateb)
    qbinvqa = qb \ qa

    BFb = szeros(T, 3, 3)
    Bτb = rotation_matrix(inv(qb)) * rotation_matrix(qa)
    return [BFb; Bτb]
end

@inline function input_jacobian_configuration_parent(joint::Rotational{T}, statea::State, stateb::State) where T
    _, qa = current_configuration(statea)
    _, qb = current_configuration(stateb)
    τ = joint.Fτ

    FaXa = szeros(T,3,3)
    FaQa = szeros(T,3,4)
    τaXa = szeros(T,3,3)
    τaQa = szeros(T,3,4)
    FbXa = szeros(T,3,3)
    FbQa = szeros(T,3,4)
    τbXa = szeros(T,3,3)
    τbQa = rotation_matrix(inv(qb)) * ∂qrotation_matrix(qa, τ)#*LVᵀmat(qa)
    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end

@inline function input_jacobian_configuration_child(joint::Rotational{T}, statea::State, stateb::State) where T
    _, qa = current_configuration(statea)
    _, qb = current_configuration(stateb)
    τ = joint.Fτ

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,4)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,4)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,4)
    τbXb = szeros(T,3,3)
    τbQb = ∂qrotation_matrix_inv(qb, rotation_matrix(qa)*τ)#*LVᵀmat(qb)
    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
