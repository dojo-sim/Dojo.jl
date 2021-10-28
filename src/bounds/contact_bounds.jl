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
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}

    n = length(normal)
    @assert n == length(bodies) == length(normal) == length(cf) == length(p) == length(offset)
    contineqcs = Vector{InequalityConstraint}()
    for i = 1:n
        contineqc = contactconstraint(bodies[i], normal[i], cf[i], p = p[i], offset = offset[i])
        push!(contineqcs, contineqc)
    end
    contineqcs = [contineqcs...] # vector typing
    return contineqcs
end

function contactconstraint(body::Body{T}, normal::AbstractVector{<:AbstractVector}, cf::AbstractVector{T};
        p::AbstractVector = [szeros(T, 3) for i=1:length(normal)],
        offset::AbstractVector = [szeros(T, 3) for i=1:length(normal)]) where {T}
    n = length(normal)
    @assert n == length(normal) == length(cf) == length(p) == length(offset)
    return contactconstraint(fill(body, n), normal, cf, p = p, offset = offset)
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
        offset::AbstractVector{T} = szeros(T, 3)) where {T}

    contbound = ContactBound(body, normal, cf, p = p)
    contineqcs = InequalityConstraint((contbound, body.id, nothing))
    return contineqcs
end

function gs(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    return g(mechanism, ineqc)
end

function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}}}
    cont = ineqc.constraints[1]
    body = getbody(mechanism, ineqc.parentid)
    x, v, q, ω = fullargssol(body.state)
    x3, q3 = posargsnext(body.state, mechanism.Δt)

    # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
    Bxmat = cont.Bx
    Bqmat = Bxmat * ∂vrotate∂q(cont.p, q3) * LVᵀmat(q3)
    SVector{4,T}(
        cont.ainv3 * (x3 + vrotate(cont.p,q3) - cont.offset) - ineqc.ssol[2][1],
        cont.cf * ineqc.γsol[2][1] - ineqc.γsol[2][2],
        (Bxmat * v + Bqmat * ω - ineqc.ssol[2][3:4])...)
end


## Derivatives accounting for quaternion specialness
## maps contact forces into the dynamics
@inline function ∂g∂pos(cont::ContactBound, x::AbstractVector, q::UnitQuaternion)
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
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    return [ineqc.γsol[2][1] * ineqc.ssol[2][1]; cone_product(ineqc.γsol[2][2:4], ineqc.ssol[2][2:4])]
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    return [ineqc.γsol[2][1] * ineqc.ssol[2][1]; cone_product(ineqc.γsol[2][2:4], ineqc.ssol[2][2:4])] - mechanism.μ * neutral_vector(ineqc.constraints[1])
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

function neutral_vector(bound::ContactBound{T,N}) where {T,N}
    N½ = Int(N/2)
    return [sones(T, 2); szeros(T, N½-2)]
end

∂g∂ʳpos(bound::ContactBound, state::State, Δt) = ∂g∂ʳpos(bound, posargsnext(state, Δt)...)

@inline function ∂g∂ʳpos(bound::ContactBound, x::AbstractVector, q::UnitQuaternion)
    X, Q = ∂g∂pos(bound, x, q)
    Q = Q # we account for quaternion specialness in ∂g∂pos
    return [X Q]
end

@inline function ∂g∂ʳvel(cont::ContactBound, x3::AbstractVector, q3::UnitQuaternion,
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

    V = X
    Ω = Q
    return [V Ω]
end

@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
    # ∇ssol[γsol .* ssol - μ; g - s] = [diag(γsol); -diag(0,1,1)]
    # ∇γsol[γsol .* ssol - μ; g - s] = [diag(ssol); -diag(1,0,0)]
    # (cf γ - ψ) dependent of ψ = γsol[2][1:1]
    # B(z) * zdot - sβ dependent of sβ = ssol[2][2:end]
    cf = ineqc.constraints[1].cf
    ∇s = [ineqc.γsol[2][1] zeros(1,3); zeros(3,1) ∇cone_product(ineqc.γsol[2][2:4]); Diagonal([-1, 0, -1, -1])]
    ∇γ = [ineqc.ssol[2][1] zeros(1,3); zeros(3,1) ∇cone_product(ineqc.ssol[2][2:4]); zeros(1,4); cf -1 0 0; zeros(2,4)]
    # matrix_entry.value = [∇s ∇γ]
    matrix_entry.value = [[ineqc.γsol[2][1] zeros(1,3); zeros(3,1) ∇cone_product(ineqc.γsol[2][2:4]); Diagonal([-1, 0, -1, -1])] [ineqc.ssol[2][1] zeros(1,3); zeros(3,1) ∇cone_product(ineqc.ssol[2][2:4]); zeros(1,4); cf -1 0 0; zeros(2,4)]]
    
    # [-γsol .* ssol + μ; -g + s]
    vector_entry.value = [-complementarityμ(mechanism, ineqc);-gs(mechanism, ineqc)]
    return
end
