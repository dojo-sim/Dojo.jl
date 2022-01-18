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
    bound = ineqc.constraints[1]
    body = getbody(mechanism, ineqc.parentid)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, mechanism.Δt)

    # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
    # vp = V(cp, B / W)_w velocity of the contact point cp, attached to body B wrt world frame, expressed in the world frame.
    vp = v25 + skew(vrotate(ϕ25, q3)) * (vrotate(bound.p, q3) - bound.offset)
    γ = ineqc.γsol[2][1]
    sγ = ineqc.ssol[2][1]
    ψ = ineqc.γsol[2][2]
    sψ = ineqc.ssol[2][2]
    β = ineqc.γsol[2][@SVector [3,4,5,6]]
    sβ = ineqc.ssol[2][@SVector [3,4,5,6]]
    SVector{6,T}(
        bound.ainv3 * (x3 + vrotate(bound.p,q3) - bound.offset) - sγ,
        bound.cf * γ - sum(β) - sψ,
        (bound.Bx * vp + ψ * sones(4) - sβ)...)
end

@inline function forcemapping(bound::LinearContactBound)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    return X
end

## Complementarity
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}; scaling = false) where {T,N,Nc,Cs<:Tuple{LinearContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    return γ .* s
end

function neutral_vector(bound::LinearContactBound{T,N}) where {T,N}
    N½ = Int(N/2)
    return sones(T, N½)
end

@inline function G(bound::LinearContactBound, x::AbstractVector, q::UnitQuaternion, λ)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    # q * ... is a rotation by quatrnon q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return [X Q]
end

@inline function ∂g∂v(bound::LinearContactBound, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, Δt
    )
    V = [bound.ainv3 * Δt;
         szeros(1,3);
         bound.Bx]
    # Ω = FiniteDiff.finite_difference_jacobian(ϕ25 -> g(bound, s, γ, x2+Δt*v25, getq3(q2,ϕ25,Δt), v25, ϕ25), ϕ25)
    ∂v∂q3 = skew(vrotate(ϕ25, q3)) * ∂vrotate∂q(bound.p, q3)
    ∂v∂q3 += skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂q(ϕ25, q3)
    ∂v∂ϕ25 = skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂p(ϕ25, q3)
    Ω = [bound.ainv3 * ∂vrotate∂q(bound.p, q3) * ∂integrator∂ϕ(q2, ϕ25, Δt)
        szeros(1,3);
        bound.Bx * (∂v∂ϕ25 + ∂v∂q3 * ∂integrator∂ϕ(q2, ϕ25, Δt))]
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

function cone_degree(bound::LinearContactBound{T,N}) where {T,N}
    # http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
    # section 2
    return Int(N / 2)
end
