mutable struct ContactBound{T,N} <: Bound{T,N}
    cf::T
    Bx::SMatrix{2,3,T,6}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}

    function ContactBound(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        Bx = SA{T}[
            1 0 0
            0 1 0
        ]

        new{Float64,8}(cf, Bx, ainv3, p, offset)
    end
end


"""
    Generate contact inequality constraints attached to a list of bodies. You need to provide:
    - the normal for each contact point
    - the coefficient of friction for each contact point
    - the offset vector p with respect to the center of the body for each contact point (optional)
    - the altitude offset for each contact point (optional)
"""
function contactconstraint(bodies::AbstractVector{<:Body{T}}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        names::AbstractVector{String} = fill("", length(normal)) ) where {T}

    n = length(normal)
    @assert n == length(bodies) == length(normal) == length(cf) == length(p) == length(offset)
    contineqcs = Vector{InequalityConstraint}()
    for i = 1:n
        contineqc = contactconstraint(bodies[i], normal[i], cf[i], p = p[i], offset = offset[i], name = names[i])
        push!(contineqcs, contineqc)
    end
    contineqcs = [contineqcs...] # vector typing
    return contineqcs
end

function contactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        names::AbstractVector{String} = fill("", length(normal)) ) where {T}
    n = length(normal)
    @assert n == length(normal) == length(cf) == length(p) == length(offset)
    return contactconstraint(fill(body, n), normal, cf, p = p, offset = offset, names = names)
end

"""
    Generate contact inequality constraint attached to one body. You need to provide:
    - the normal for the contact point
    - the coefficient of friction for the contact point
    - the offset vector p with respect to the center of the body for the contact point (optional)
    - the altitude offset for the contact point (optional)
"""
function contactconstraint(body::Body{T}, normal::AbstractVector{T}, cf::T;
        p::AbstractVector{T} = szeros(T, 3),
        offset::AbstractVector{T} = szeros(T, 3), name::String = "") where {T}

    contbound = ContactBound(body, normal, cf, p = p, offset = offset)
    contineqcs = InequalityConstraint((contbound, body.id, nothing); name = name)
    return contineqcs
end

function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}}}
    bound = ineqc.constraints[1]
    body = getbody(mechanism, ineqc.parentid)
    x2, v25, q2, ϕ25 = fullargssol(body.state)
    x3, q3 = posargs3(body.state, mechanism.Δt)

    g(bound, ineqc.ssol[2], ineqc.γsol[2], x3, q3, v25, ϕ25)
end

function g(bound::ContactBound, s::AbstractVector{T}, γ::AbstractVector{T},
        x3::AbstractVector{T}, q3::UnitQuaternion{T}, v25::AbstractVector{T},
        ϕ25::AbstractVector{T}) where {T}

    # transforms the velocities of the origin of the link into velocities
    vp = v25 + skew(vrotate(ϕ25, q3)) * (vrotate(bound.p, q3) - bound.offset)
    SVector{4,T}(
        bound.ainv3 * (x3 + vrotate(bound.p, q3) - bound.offset) - s[1],
        bound.cf * γ[1] - γ[2],
        (bound.Bx * vp - s[@SVector [3,4]])...)
end


## Derivatives accounting for quaternion specialness
## maps contact forces into the dynamics
@inline function G(bound::ContactBound, x::AbstractVector, q::UnitQuaternion, λ)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    # q * ... is a rotation by quatrnon q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return [X Q]
end

@inline function forcemapping(bound::ContactBound)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    return X
end

## Complementarity
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}; scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    if scaling
        W, Wi, λ = nt_scaling(ineqc)
        c = vcat(λ[1] * λ[1], cone_product(λ[@SVector [2,3,4]], λ[@SVector [2,3,4]]))
    else
        c = vcat(γ[1] * s[1], cone_product(γ[@SVector [2,3,4]], s[@SVector [2,3,4]]))
    end
    return c
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}; scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    if scaling
        W, Wi, λ = nt_scaling(ineqc)
        c = vcat(λ[1] * λ[1], cone_product(λ[@SVector [2,3,4]], λ[@SVector [2,3,4]])) - mechanism.μ * neutral_vector(ineqc.constraints[1])
    else
        c = vcat(γ[1] * s[1], cone_product(γ[@SVector [2,3,4]], s[@SVector [2,3,4]])) - mechanism.μ * neutral_vector(ineqc.constraints[1])
    end
    return c
end

function neutral_vector(bound::ContactBound{T,N}) where {T,N}
    N½ = Int(N/2)
    return [sones(T, 2); szeros(T, N½-2)]
end

@inline function ∂g∂ʳvel(bound::ContactBound{T}, x3::AbstractVector{T}, q3::UnitQuaternion{T},
        x2::AbstractVector{T}, v25::AbstractVector{T}, q2::UnitQuaternion{T}, ϕ25::AbstractVector{T}, λ, Δt::T) where {T}
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


@inline function ∂g∂ʳvel(bound::ContactBound, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, Δt
    )
    Bxmat = bound.Bx
    p = bound.p

    V = [bound.ainv3 * Δt;
         szeros(1,3);
         Bxmat]

    ∂v∂q3 = ∂vrotate∂q(-skew(p) * ϕ25, q3)
    ∂v∂ϕ25 = ∂vrotate∂p(-skew(p) * ϕ25, q3) * -skew(p)
    Ω = [bound.ainv3 * ∂vrotate∂q(p, q3) * ∂integrator∂ϕ(q2, ϕ25, Δt)
        szeros(1,3);
        Bxmat * (∂v∂ϕ25 + ∂v∂q3 * ∂integrator∂ϕ(q2, ϕ25, Δt))]
    return [V Ω]
end

@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    # ∇ssol[γsol .* ssol - μ; g - s] = [diag(γsol); -diag(0,1,1)]
    # ∇γsol[γsol .* ssol - μ; g - s] = [diag(ssol); -diag(1,0,0)]
    # (cf γ - ψ) dependent of ψ = γsol[2][1:1]
    # B(z) * zdot - sβ dependent of sβ = ssol[2][2:end]
    cf = ineqc.constraints[1].cf
    γ = ineqc.γsol[2] + 1e-10*neutral_vector(ineqc.constraints[1]) # TODO need to check this is legit
    s = ineqc.ssol[2] + 1e-10*neutral_vector(ineqc.constraints[1]) # TODO need to check this is legit

    # ∇s = [ineqc.γsol[2][1] szeros(1,3); szeros(3,1) ∇cone_product(ineqc.γsol[2][2:4]); Diagonal([-1, 0, -1, -1])]
    ∇s1 = [γ[SA[1]]; szeros(T,3)]'
    ∇s2 = [szeros(T,3,1) ∇cone_product(γ[SA[2,3,4]])]
    ∇s3 = Diagonal(SVector{4,T}(-1, 0, -1, -1))
    ∇s = [∇s1; ∇s2; ∇s3]

    # ∇γ = [ineqc.ssol[2][1] szeros(1,3); szeros(3,1) ∇cone_product(ineqc.ssol[2][2:4]); szeros(1,4); cf -1 0 0; szeros(2,4)]
    ∇γ1 = [s[SA[1]]; szeros(T,3)]'
    ∇γ2 = [szeros(T,3,1) ∇cone_product(s[SA[2,3,4]])]
    ∇γ3 = SA[0   0 0 0;
             cf -1 0 0;
             0   0 0 0;
             0   0 0 0;]
    ∇γ = [∇γ1; ∇γ2; ∇γ3]

    # matrix_entry.value = [[ineqc.γsol[2][1] szeros(1,3); szeros(3,1) ∇cone_product(ineqc.γsol[2][2:4]); Diagonal([-1, 0, -1, -1])] [ineqc.ssol[2][1] szeros(1,3); szeros(3,1) ∇cone_product(ineqc.ssol[2][2:4]); szeros(1,4); cf -1 0 0; szeros(2,4)]]
    matrix_entry.value = [∇s ∇γ]

    # [-γsol .* ssol + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, ineqc), -g(mechanism, ineqc))
    return
end

function nt_scaling(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    W_ort, Wi_ort, λ_ort = ort_nt_scaling(s[1:1], γ[1:1])
    W_soc, Wi_soc, λ_soc = soc_nt_scaling(s[2:4], γ[2:4])
    W = cat(W_ort, W_soc, dims = (1,2))
    Wi = cat(Wi_ort, Wi_soc, dims = (1,2))
    λ = [λ_ort; λ_soc]
    return W, Wi, λ
end

function cone_degree(bound::ContactBound{T,N}) where {T,N}
    # http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
    # section 2
    return 1 + 1 # 1 for impact and 1 for the second order cone
end
