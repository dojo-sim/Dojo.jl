mutable struct ImpactBound{T,N} <: Bound{T,N}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}

    function ImpactBound(body::Body{T}, normal::AbstractVector; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(ainv3, p, offset)
    end
end

function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ImpactBound{T,N}}}
    bound = ineqc.constraints[1]
    body = getbody(mechanism, ineqc.parentid)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, mechanism.Δt)
    SVector{1,T}(bound.ainv3 * (x3 + vrotate(bound.p,q3) - bound.offset) - ineqc.ssol[2][1])
end

@inline function ∂g∂v(bound::ImpactBound, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, Δt)
    V = bound.ainv3 * Δt
    Ω = bound.ainv3 * ∂vrotate∂q(bound.p, q3) * ∂integrator∂ϕ(q2, ϕ25, Δt)
    return [V Ω]
end

@inline function G(bound::ImpactBound, x::AbstractVector, q::UnitQuaternion, λ)
    X = bound.ainv3
    # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return [X Q]
end

@inline function forcemapping(bound::ImpactBound)
    X = bound.ainv3
    return X
end

@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ImpactBound{T,N}},N½}
    # ∇ssol[γsol .* ssol - μ; g - s] = [diag(γsol); -diag(0,1,1)]
    # ∇γsol[γsol .* ssol - μ; g - s] = [diag(ssol); -diag(1,0,0)]
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-γsol .* ssol + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, ineqc), -g(mechanism, ineqc))
    return
end

"""
    Generate impact inequality constraints attached to a list of bodies. You need to provide:
    - the normal for each contact point
    - the coefficient of friction for each contact point
    - the offset vector p with respect to the center of the body for each contact point (optional)
    - the altitude offset for each contact point (optional)
"""
function impactconstraint(bodies::AbstractVector{<:Body{T}}, normal::AbstractVector{<:AbstractVector};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}

    n = length(normal)
    @assert n == length(bodies) == length(normal) == length(p) == length(offset)
    impineqcs = Vector{InequalityConstraint}()
    for i = 1:n
        impineqc = impactconstraint(bodies[i], normal[i], p = p[i], offset = offset[i])
        push!(impineqcs, impineqc)
    end
    impineqcs = [impineqcs...] # vector typing
    return impineqcs
end

function impactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
    n = length(normal)
    @assert n == length(normal) == length(p) == length(offset)
    return impactconstraint(fill(body, n), normal, p = p, offset = offset)
end

"""
    Generate impact inequality constraint attached to one body. You need to provide:
    - the normal for the contact point
    - the offset vector p with respect to the center of the body for the contact point (optional)
    - the altitude offset for the contact point (optional)
"""
function impactconstraint(body::Body{T}, normal::AbstractVector{T};
        p::AbstractVector{T} = szeros(T, 3),
        offset::AbstractVector{T} = szeros(T, 3)) where {T}

    impbound = ImpactBound(body, normal, p = p, offset = offset)
    impineqcs = InequalityConstraint((impbound, body.id, nothing))
    return impineqcs
end