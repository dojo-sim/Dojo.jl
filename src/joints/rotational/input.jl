@inline function applyFτ!(joint::Rotational{T}, statea::State, stateb::State, Δt::T, clear::Bool) where T
    τ = joint.Fτ
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)

    τa = vrotate(-τ, qa) # in world coordinates
    τb = -τa # in world coordinates

    τa = vrotate(τa,inv(qa)) # in local coordinates
    τb = vrotate(τb,inv(qb)) # in local coordinates

    statea.τ2[end] += τa
    stateb.τ2[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

@inline function ∂Fτ∂ua(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    BFa = (szeros(T, 3, 3))
    Bτa = -I

    return [BFa; Bτa]
end

@inline function ∂Fτ∂ub(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)
    qbinvqa = qb \ qa

    BFb = (szeros(T, 3, 3))
    Bτb = VLmat(qbinvqa) * RᵀVᵀmat(qbinvqa)

    return [BFb; Bτb]
end

@inline function ∂Fτ∂a(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)
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

@inline function ∂Fτ∂b(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)
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
