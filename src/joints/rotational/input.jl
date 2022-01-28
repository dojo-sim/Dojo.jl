@inline function apply_input!(joint::Rotational{T}, statea::State, stateb::State, timestep::T, clear::Bool) where T
    τ = joint.Fτ
    _, qa = current_configuration(statea)
    _, qb = current_configuration(stateb)

    τa = vrotate(-τ, qa) # in world coordinates
    τb = -τa # in world coordinates

    τa = vrotate(τa,inv(qa)) # in local coordinates
    τb = vrotate(τb,inv(qb)) # in local coordinates

    statea.τ2[end] += τa
    stateb.τ2[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

@inline function input_jacobian_control_parent(joint::Rotational{T}, statea::State, stateb::State, timestep::T) where T
    BFa = (szeros(T, 3, 3))
    Bτa = -I

    return [BFa; Bτa]
end

@inline function input_jacobian_control_child(joint::Rotational{T}, statea::State, stateb::State, timestep::T) where T
    _, qa = current_configuration(statea)
    _, qb = current_configuration(stateb)
    qbinvqa = qb \ qa

    BFb = (szeros(T, 3, 3))
    Bτb = VLmat(qbinvqa) * RᵀVᵀmat(qbinvqa)

    return [BFb; Bτb]
end

@inline function input_jacobian_configuration_parent(joint::Rotational{T}, statea::State, stateb::State, timestep::T) where T
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
    τbQa = 2*VLᵀmat(qb)*Rmat(qb)*Rᵀmat(qa)*Rmat(UnitQuaternion(τ))#*LVᵀmat(qa)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end

@inline function input_jacobian_configuration_child(joint::Rotational{T}, statea::State, stateb::State, timestep::T) where T
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
    τbQb = 2*VLᵀmat(qb)*Lmat(qa)*Lmat(UnitQuaternion(τ))*Lᵀmat(qa)#*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end