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

@inline function ∂g∂v(bound::ContactBound{T}, x3::AbstractVector{T}, q3::UnitQuaternion{T},
    x2::AbstractVector{T}, v25::AbstractVector{T}, q2::UnitQuaternion{T}, ϕ25::AbstractVector{T}, λ, Δt::T) where {T}
    V = [bound.ainv3 * Δt;
        szeros(1,3);
        bound.Bx]
    # Ω = FiniteDiff.finite_difference_jacobian(ϕ25 -> g(bound, s, γ, x2+Δt*v25, getq3(q2,ϕ25,Δt), v25, ϕ25), ϕ25)
    ∂v∂q3 = skew(vrotate(ϕ25, q3)) * ∂vrotate∂q(bound.p, q3)
    ∂v∂q3 += skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂q(ϕ25, q3)
    ∂v∂ϕ25 = skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂p(ϕ25, q3)
    Ω = [bound.ainv3 * ∂vrotate∂q(bound.p, q3) * ∂integrator∂ϕ(q2, ϕ25, Δt);
        szeros(1,3);
        bound.Bx * (∂v∂ϕ25 + ∂v∂q3 * ∂integrator∂ϕ(q2, ϕ25, Δt))]
    return [V Ω]
end

@inline function ∂g∂z(bound::ContactBound{T}, x3::AbstractVector{T}, q3::UnitQuaternion{T},
    x2::AbstractVector{T}, v25::AbstractVector{T}, q2::UnitQuaternion{T}, ϕ25::AbstractVector{T}, λ, Δt::T) where {T}
    X = [bound.ainv3;
        szeros(1,3);
        szeros(2,3)]
    # Ω = FiniteDiff.finite_difference_jacobian(ϕ25 -> g(bound, s, γ, x2+Δt*v25, getq3(q2,ϕ25,Δt), v25, ϕ25), ϕ25)
    ∂v∂q3 = skew(vrotate(ϕ25, q3)) * ∂vrotate∂q(bound.p, q3)
    ∂v∂q3 += skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂q(ϕ25, q3)
    Q = [bound.ainv3 * ∂vrotate∂q(bound.p, q3);
        szeros(1,4);
        bound.Bx * ∂v∂q3]
    return [X Q]
end

@inline function G(bound::ContactBound, x::AbstractVector, q::UnitQuaternion, λ)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return [X Q]
end

@inline function forcemapping(bound::ContactBound)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    return X
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

function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½};
        scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    γ = ineqc.γsol[2]
    s = ineqc.ssol[2]
    return vcat(γ[1] * s[1], cone_product(γ[@SVector [2,3,4]], s[@SVector [2,3,4]]))
end

neutral_vector(bound::ContactBound{T,N}) where {T,N} = [sones(T, 2); szeros(T, Int(N/2) -2)]

cone_degree(bound::ContactBound) = 2

## Utilities
# signed distance function
function sdf(ineqc::InequalityConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
        q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    cont = ineqc.constraints[1]
    return cont.ainv3 * (x + vrotate(cont.p, q) - cont.offset)
end

function get_sdf(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, storage::Storage{T,N}) where {T,Nn,Ne,Nb,Ni,N}
    d = []
    for ineqc in mechanism.ineqconstraints
        ibody = getbody(mechanism, ineqc.parentid).id - Ne
        push!(d, [sdf(ineqc, storage.x[ibody][i], storage.q[ibody][i]) for i = 1:N])
    end
    return d
end

# contact location
function contact_location(mechanism::Mechanism)
    return [contact_location(mech, ineqc) for ineqc in mechanism.ineqconstraints]
end

function contact_location(mechanism::Mechanism, ineqc::InequalityConstraint)
    body = mechanism.bodies[findfirst(x -> x.id == ineqc.parentid, mechanism.bodies)]
    return contact_location(ineqc, body)
end

function contact_location(ineqc::InequalityConstraint{T,N,Nc,Cs},
        body::Body) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    x = body.state.x2[1]
    q = body.state.q2[1]
    return contact_location(ineqc, x, q)
end

function contact_location(ineqc::InequalityConstraint{T,N,Nc,Cs}, x::AbstractVector{T},
        q::UnitQuaternion{T}) where {T,N,Nc,Cs<:Tuple{<:Bound{T,N}}}
    cont = ineqc.constraints[1]
    return x + vrotate(cont.p,q) - cont.offset
end
