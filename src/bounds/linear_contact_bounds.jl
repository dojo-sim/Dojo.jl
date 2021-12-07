mutable struct LinearContactBound{T,N} <: Bound{T,N}
    cf::T
    Bx::SMatrix{4,3,T,12}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}

    function LinearContactBound(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        Bx = SA{T}[
             1  0  0
            -1  0  0
             0  1  0
             0 -1  0
        ]
        new{Float64,12}(cf, Bx, ainv3, p, offset)
    end
end

"""
    Generate linear contact inequality constraints attached to a list of bodies. You need to provide:
    - the normal for each contact point
    - the coefficient of friction for each contact point
    - the offset vector p with respect to the center of the body for each contact point (optional)
    - the altitude offset for each contact point (optional)
"""
function linearcontactconstraint(bodies::AbstractVector{<:Body{T}}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}

    n = length(normal)
    @assert n == length(bodies) == length(normal) == length(cf) == length(p) == length(offset)
    linineqcs = Vector{InequalityConstraint}()
    for i = 1:n
        linineqc = linearcontactconstraint(bodies[i], normal[i], cf[i], p = p[i], offset = offset[i])
        push!(linineqcs, linineqc)
    end
    linineqcs = [linineqcs...] # vector typing
    return linineqcs
end

function linearcontactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
    n = length(normal)
    @assert n == length(normal) == length(cf) == length(p) == length(offset)
    return linearcontactconstraint(fill(body, n), normal, cf, p = p, offset = offset)
end

"""
    Generate linear contact inequality constraint attached to one body. You need to provide:
    - the normal for the contact point
    - the coefficient of friction for the contact point
    - the offset vector p with respect to the center of the body for the contact point (optional)
    - the altitude offset for the contact point (optional)
"""
function linearcontactconstraint(body::Body{T}, normal::AbstractVector{T}, cf::T;
        p::AbstractVector{T} = szeros(T, 3),
        offset::AbstractVector{T} = szeros(T, 3)) where {T}

    linbound = LinearContactBound(body, normal, cf, p = p, offset = offset)
    linineqcs = InequalityConstraint((linbound, body.id, nothing))
    return linineqcs
end

function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{LinearContactBound{T,N}}}
    cont = ineqc.constraints[1]
    body = getbody(mechanism, ineqc.parentid)
    x, v, q, ω = fullargssol(body.state)
    x3, q3 = posargs3(body.state, mechanism.Δt)

    # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
    Bxmat = cont.Bx
    Bqmat = Bxmat * ∂vrotate∂q(cont.p, q3) * LVᵀmat(q3)
    γ = ineqc.γsol[2][1]
    sγ = ineqc.ssol[2][1]
    ψ = ineqc.γsol[2][2]
    sψ = ineqc.ssol[2][2]
    β = ineqc.γsol[2][@SVector [3,4,5,6]]
    sβ = ineqc.ssol[2][@SVector [3,4,5,6]]
    SVector{6,T}(
        cont.ainv3 * (x3 + vrotate(cont.p,q3) - cont.offset) - sγ,
        cont.cf * γ - sum(β) - sψ,
        (Bxmat * v + Bqmat * ω + ψ * sones(4) - sβ)...)
end


## Derivatives accounting for quaternion specialness
## maps contact forces into the dynamics
@inline function ∂g∂pos(cont::LinearContactBound, x::AbstractVector, q::UnitQuaternion)
    Bxmat = cont.Bx
    p = cont.p
    nx = size(x)[1]
    nq = nx

    drot = ∂vrotate∂q(cont.p, q) * LVᵀmat(q)

    X = [cont.ainv3;
         szeros(1,nx);
         Bxmat]
    Q = [cont.ainv3 * drot;
         szeros(1,nq);
         Bxmat * drot]
    return X, Q
end


## Complementarity
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}; scaling = false) where {T,N,Nc,Cs<:Tuple{LinearContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    if scaling
        W, Wi, λ = nt_scaling(ineqc)
        c = λ .* λ
    else
        c = γ .* s
    end
    return c
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}; scaling = false) where {T,N,Nc,Cs<:Tuple{LinearContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    if scaling
        W, Wi, λ = nt_scaling(ineqc)
        c = λ .* λ - mechanism.μ * neutral_vector(ineqc.constraints[1])
    else
        c = γ .* s - mechanism.μ * neutral_vector(ineqc.constraints[1])
    end
    return c
end

function neutral_vector(bound::LinearContactBound{T,N}) where {T,N}
    N½ = Int(N/2)
    return sones(T, N½)
end

@inline function ∂g∂ʳpos(bound::LinearContactBound, x::AbstractVector, q::UnitQuaternion)
    X, Q = ∂g∂pos(bound, x, q)
    Q = Q # we account for quaternion specialness in ∂g∂pos
    return [X Q]
end

@inline function ∂g∂ʳvel(cont::LinearContactBound, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v2::AbstractVector, q2::UnitQuaternion, ω2::AbstractVector, Δt
    )
    Bxmat = cont.Bx
    p = cont.p
    nx = size(x2)[1]
    nq = nx


    X = [cont.ainv3 * Δt;
         szeros(1,nx);
         Bxmat]

    B(q) = Bxmat * ∂vrotate∂q(cont.p, UnitQuaternion(q...)) * LVᵀmat(UnitQuaternion(q...))

    Q = [(cont.ainv3 * (VLmat(q3) * Lmat(UnitQuaternion(cont.p)) * Tmat() + VRᵀmat(q3) * Rmat(UnitQuaternion(cont.p)))) * Lmat(q2) * derivωbar(ω2, Δt) * Δt / 2
         szeros(1,nq);
         B(q3) + ForwardDiff.jacobian(q -> B(q) * ω2, [q3.w; q3.x; q3.y; q3.z]) * Lmat(q2) * derivωbar(ω2, Δt) * Δt / 2]
         # B(q3) + FiniteDiff.finite_difference_jacobian(q -> B(q) * ω2, [q3.w; q3.x; q3.y; q3.z]) * Lmat(q2) * derivωbar(ω2, Δt) * Δt / 2]

    V = X
    Ω = Q
    return [V Ω]
end

@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{LinearContactBound{T,N}},N½}
    # ∇ssol[γsol .* ssol - μ; g - s] = [diag(γsol); -diag(0,1,1)]
    # ∇γsol[γsol .* ssol - μ; g - s] = [diag(ssol); -diag(1,0,0)]
    # (cf γ - ψ) dependent of ψ = γsol[2][1:1]
    # B(z) * zdot - sβ dependent of sβ = ssol[2][2:end]
    cf = ineqc.constraints[1].cf
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]

    ∇s1 = Diagonal(γ) # 6x6
    ∇s2 = Diagonal(-sones(T,6))
    ∇s = vcat(∇s1, ∇s2) # 12x6

    ∇γ1 = Diagonal(s) # 6x6
    ∇γ2 = @SMatrix[ 0  0  0  0  0  0;
                   cf  0 -1 -1 -1 -1;
                    0  1  0  0  0  0;
                    0  1  0  0  0  0;
                    0  1  0  0  0  0;
                    0  1  0  0  0  0;]
    ∇γ = vcat(∇γ1, ∇γ2) # 12x6
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-γsol .* ssol + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, ineqc), -g(mechanism, ineqc))
    return
end

function nt_scaling(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{LinearContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    return ort_nt_scaling(s, γ)
end

function cone_degree(bound::LinearContactBound{T,N}) where {T,N}
    # http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
    # section 2
    return Int(N/2)
end
