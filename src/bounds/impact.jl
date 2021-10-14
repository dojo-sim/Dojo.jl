mutable struct Impact{T,N} <: Bound{T,N}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}


    function Impact(body::Body{T}, normal::AbstractVector; p::AbstractVector = zeros(3),  offset::AbstractVector = zeros(3)) where T
        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        # p is the 3d location of the contact point attached to the link in the link's frame
        # offset is the 3d location of the point on the surface closest to the contact point
        # v3 is the normal of
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'

        new{T,2}(ainv3, p, offset), body.id, nothing
    end
end

### Constraints and derivatives
## Position level constraints (for dynamics)
function g(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{Impact{T,N}}}
    g(ineqc.constraints[1], getbody(mechanism, ineqc.parentid), mechanism.Δt)
end

@inline g(impact::Impact{T}, x::AbstractVector, q::UnitQuaternion) where T = SVector{1,T}(impact.ainv3 * (x + vrotate(impact.p,q) - impact.offset))

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂pos(impact::Impact, x::AbstractVector, q::UnitQuaternion)
    p = impact.p
    X = impact.ainv3
    Q = impact.ainv3 * (VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p)))
    return X, Q
end

## Complementarity
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{Impact{T,N}},N½}
    return ineqc.γsol[2] .* ineqc.ssol[2]
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{Impact{T,N}},N½}
    return ineqc.γsol[2] .* ineqc.ssol[2] .- mechanism.μ
end

function neutral_vector(bound::Bound{T,N}) where {T,N}
    N½ = Int(N/2)
    return sones(T, N½)
end