mutable struct ConeBound{T,N} <: Bound{T,N}
    cf::T
    Bx::SMatrix{2,3,T,6}
    p::SVector{3,T}
    γsolref::Vector{SVector{1,T}} #TODO this is a reference to the associated impact γ to avoid allocations

    function ConeBound(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T

        Bx = SA{T}[
            1 0 0
            0 1 0
        ]

        impact = InequalityConstraint(Impact(body, normal; p = p, offset = offset))
        γsolref = impact.γsol
        new{Float64,6}(cf, Bx, p, γsolref), impact
    end
end

"""
    Generate cone and impact inequality constraints attached to a list of bodies. You need to provide:
    - the normal for each contact point
    - the coefficient of friction for each contact point
    - the offset vector p with respect to the center of the body for each contact point (optional)
    - the altitude offset for each contact point (optional)
"""
function splitcontactconstraint(bodies::AbstractVector{<:Body{T}}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}

    n = length(normal)
    @assert n == length(bodies) == length(normal) == length(cf) == length(p) == length(offset)
    conineqcs = Vector{InequalityConstraint}()
    impineqcs = Vector{InequalityConstraint}()
    for i = 1:n
        conineqc, impineqc = splitcontactconstraint(bodies[i], normal[i], cf[i], p = p[i], offset = offset[i])
        push!(conineqcs, conineqc)
        push!(impineqcs, impineqc)
    end
    conineqcs = [conineqcs...] # vector typing
    impineqcs = [impineqcs...] # vector typing
    return conineqcs, impineqcs
end

function splitcontactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
    n = length(normal)
    @assert n == length(normal) == length(cf) == length(p) == length(offset)
    return splitcontactconstraint(fill(body, n), normal, cf, p = p, offset = offset)
end

"""
    Generate cone and impact inequality constraint attached to one body. You need to provide:
    - the normal for the contact point
    - the coefficient of friction for the contact point
    - the offset vector p with respect to the center of the body for the contact point (optional)
    - the altitude offset for the contact point (optional)
"""
function splitcontactconstraint(body::Body{T}, normal::AbstractVector{T}, cf::T;
        p::AbstractVector{T} = szeros(T, 3),
        offset::AbstractVector{T} = szeros(T, 3)) where {T}

    tmp = ConeBound(body, normal, cf; p = p)
    conbound = getindex(tmp,1)
    impineqc = getindex(tmp,2)
    impid = getfield(impineqc, :id)
    conineqc = InequalityConstraint((conbound, body.id, impid))
    return conineqc, impineqc
end

function gs(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ConeBound{T,N}},N½}
    # this is the residual with substracted slacks
    # error()
    # we remove the - ineqc.ssol[2] because this is not true for ConeBound
    # we already account for the - ψ and - sβ in g
    return g(mechanism, ineqc)# - ineqc.ssol[2]
end

function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ConeBound{T,N}}}
    cone = ineqc.constraints[1]
    body = getbody(mechanism, ineqc.parentid)
    x, v, q, ω = fullargssol(body.state)

    # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
    Bxmat = cone.Bx
    # transforms the velocities of the contact point attached to the link into velocities along all 4 axes of the friction pyramid
    Bqmat = Bq(Bxmat, cone.p, q)
    SVector{3,T}(
        cone.cf * cone.γsolref[2][1] - ineqc.γsol[2][1],
        (Bxmat*v + Bqmat*ω - ineqc.ssol[2][2:3])...)
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂pos(cone::ConeBound, x::AbstractVector, q::UnitQuaternion)
    Bxmat = cone.Bx
    nx = size(x)[1]
    nq = nx
    X = [szeros(1,nx); Bxmat]
    Q = [szeros(1,nq); Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-cone.p)*2.0]
    return X, Q
end

@inline function ∂gab∂ʳba(impact::Impact, cone::ConeBound{T}) where T
    SA{T}[0 0 0 0 0 0; 0 0 0 0.0 0 0], SA{T}[0 0; 0 0; 0 0; 0 cone.cf; 0 0; 0 0]
end

## Complementarity
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ConeBound{T,N}},N½}
    return cone_product(ineqc.γsol[2], ineqc.ssol[2])
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ConeBound{T,N}},N½}
    return cone_product(ineqc.γsol[2], ineqc.ssol[2]) - mechanism.μ * neutral_vector(ineqc.constraints[1])
end

function ∇cone_product(u::AbstractVector{T}) where {T}
    n = length(u)
    U = zeros(n,n)
    U += u[1] * I(n)
    U[1,2:end] += u[2:end]
    U[2:end, 1] += u[2:end]
    return U
end

function cone_product(u::AbstractVector{T}, v::AbstractVector{T}) where {T}
    [u'*v; u[1] * v[2:end] + v[1] * u[2:end]]
end

function neutral_vector(bound::ConeBound{T,N}) where {T,N}
    N½ = Int(N/2)
    return [sones(T, 1); szeros(T, N½-1)]
end

@inline function ∂g∂ʳpos(bound::ConeBound, x::AbstractVector, q::UnitQuaternion)
    X, Q = ∂g∂pos(bound, x, q)
    Q = Q# * LVᵀmat(q)
    return [X Q]
end

@inline function ∂g∂ʳvel(bound::ConeBound, x2::AbstractVector, q2::UnitQuaternion,
    x1::AbstractVector, v1::AbstractVector, q1::UnitQuaternion, ω1::AbstractVector, Δt
    )

    X, Q = ∂g∂pos(bound, x1, q1)
    V = X #* Δt
    Ω = Q #* Lmat(q1) * derivωbar(ω1, Δt) * Δt / 2

    return [V Ω]
end

@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ConeBound{T,N}},N½}
    # ∇ssol[γsol .* ssol - μ; g - s] = [diag(γsol); -diag(0,1,1)]
    # ∇γsol[γsol .* ssol - μ; g - s] = [diag(ssol); -diag(1,0,0)]
    # (cf γ - ψ) dependent of ψ = γsol[2][1:1]
    # B(z) * zdot - sβ dependent of sβ = ssol[2][2:end]
    matrix_entry.value = [[∇cone_product(ineqc.γsol[2]); zeros(1, 3); zeros(2,1) -Diagonal(ones(2))] [∇cone_product(ineqc.ssol[2]); -1.0 zeros(1, 2); zeros(2, 3)]]
    # [-γsol .* ssol + μ; -g + s]
    vector_entry.value = [-complementarityμ(mechanism, ineqc);-gs(mechanism, ineqc)]
    return
end
