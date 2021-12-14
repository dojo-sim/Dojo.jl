mutable struct Translational{T,N,N̄,Nl} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
    spring::T
    damper::T
    spring_offset::SVector{N̄,T}
    joint_limits::Vector{SVector{Nl,T}} # lower and upper limits on the joint minimal coordinate angles
    Fτ::SVector{3,T}
    function Translational{T,N,N̄}(body1::Component, body2::Component;
            p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3),
            spring = zero(T), damper = zero(T), spring_offset = szeros(T,N̄),
            joint_limits = [szeros(T,0), szeros(T,0)],
        ) where {T,N,N̄,Nl}
        vertices = (p1, p2)
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]
        Fτ = zeros(T,3)
        nl = length(joint_limits[1])
        new{T,N,N̄,nl}(V3, V12, vertices, spring, damper, spring_offset, joint_limits, Fτ), body1.id, body2.id
    end
end

joint_limits_length(joint::Translational{T,N,N̄,Nl}) where {T,N,N̄,Nl} = Nl

Translational0{T,Nl} = Translational{T,0,3,Nl} where {T,Nl}
Translational1{T,Nl} = Translational{T,1,2,Nl} where {T,Nl}
Translational2{T,Nl} = Translational{T,2,1,Nl} where {T,Nl}
Translational3{T,Nl} = Translational{T,3,0,Nl} where {T,Nl}

@inline function _g(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end

@inline function _g(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end

# Position level constraints (for dynamics)
@inline function g(joint::Translational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    vertices = joint.vertices
    e = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
    return constraintmat(joint) * e
end

@inline function g(joint::Translational{T,N,N̄,0}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    vertices = joint.vertices
    e = xb + vrotate(vertices[2], qb) - vertices[1]
    return constraintmat(joint) * e
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Translational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    point2 = xb + vrotate(joint.vertices[2], qb)
    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = ∂vrotate∂q(point2 - (xa + vrotate(joint.vertices[1], qa)), inv(qa)) * Tmat()
    Q += ∂vrotate∂p(point2 - (xa + vrotate(joint.vertices[1], qa)), inv(qa)) * -∂vrotate∂q(joint.vertices[1], qa)
    return constraintmat(joint) * X, constraintmat(joint) * Q
end

@inline function ∂g∂posb(joint::Translational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))
    return constraintmat(joint) * X, constraintmat(joint) * Q
end
@inline function ∂g∂posb(joint::Translational{T,N,N̄,0}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = I
    Q = 2 * VRᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))
    return constraintmat(joint) * X, constraintmat(joint) * Q
end

function ∂g∂ʳposa(joint::Translational{T,N,N̄,0}, statea::State, stateb::State, λ, Δt) where {T,N,N̄}
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)
    ∂g∂ʳposa(joint, xa, qa, xb, qb, λ)
end

function ∂g∂ʳposa(joint::Translational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = -1.0 * transpose(rotation_matrix(qa))
    pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    ca_a = rotation_matrix(inv(qa)) * (xa) # body a com
    capb_a = pb_a - ca_a
    Q = - 1.0 * transpose(skew(capb_a))
    return constraintmat(joint) * [X Q]
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,0}, statea::State, stateb::State, λ, Δt) where {T,N,N̄}
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)
    ∂g∂ʳposb(joint, xa, qa, xb, qb, λ)
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = transpose(rotation_matrix(qa))
    pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    cbpb_b = pb_b - cb_b
    Q = transpose(skew(cbpb_b) * rotation_matrix(inv(qb) * qa))
    return constraintmat(joint) * [X Q]
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,0}, stateb::State, λ, Δt) where {T,N,N̄}
    xb, qb = posargs2(stateb)
    ∂g∂ʳposb(joint, xb, qb, λ)
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,0}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = transpose(I(3))
    pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    cbpb_b = pb_b - cb_b
    Q = transpose(skew(cbpb_b) * rotation_matrix(inv(qb)))
    return constraintmat(joint) * [X Q]
end

### w/ Limits
function get_sγ(joint::Translational{T,N,N̄,Nl}, λ) where {T,N,N̄,Nl}
    s = λ[N .+ (1:(2 * Nl))]
    γ = λ[N + 2 * Nl .+ (1:(2 * Nl))]
    return s, γ
end

# Position level constraints (for dynamics)
@inline function g(joint::Translational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    vertices = joint.vertices
    e1 = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
    e2 = minimalCoordinates(joint, xa, qa, xb, qb)

    s, γ = get_sγ(joint, λ)

    return [
            constraintmat(joint) * e1;
            s[1:Nl] - (joint.joint_limits[2] - e2);
            s[Nl .+ (1:Nl)] - (e2 - joint.joint_limits[1]);
           ]
end

@inline function g(joint::Translational{T,N,N̄,Nl}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    vertices = joint.vertices
    e1 = xb + vrotate(vertices[2], qb) - vertices[1]
    e2 = minimalCoordinates(joint, xa, qa, xb, qb)

    s, γ = get_sγ(joint, λ)

    return [
            constraintmat(joint) * e1;
            s[1:Nl] - (joint.joint_limits[2] - e2);
            s[Nl .+ (1:Nl)] - (e2 - joint.joint_limits[1]);
           ]
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Translational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    # point2 = xb + vrotate(joint.vertices[2], qb)
    # X = -VLᵀmat(qa) * RVᵀmat(qa)
    # Q = ∂vrotate∂q(point2 - (xa + vrotate(joint.vertices[1], qa)), inv(qa)) * Tmat()
    # Q += ∂vrotate∂p(point2 - (xa + vrotate(joint.vertices[1], qa)), inv(qa)) * -∂vrotate∂q(joint.vertices[1], qa)
    X = FiniteDiff.finite_difference_jacobian(x -> g(joint, x, qa, xb, qb, λ), xa)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, UnitQuaternion(q..., false), xb, qb, λ), vector(qa))
    return X, Q
end

@inline function ∂g∂posb(joint::Translational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    # X = VLᵀmat(qa) * RVᵀmat(qa)
    # Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))
    X = FiniteDiff.finite_difference_jacobian(x -> g(joint, xa, qa, x, qb, λ), xb)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, qa, xb, UnitQuaternion(q..., false), λ), vector(qb))
    return X, Q
end
@inline function ∂g∂posb(joint::Translational{T,N,N̄,Nl}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    # X = I
    # Q = 2 * VRᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))
    X = FiniteDiff.finite_difference_jacobian(x -> g(joint, x, qb, λ), xb)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xb, UnitQuaternion(q..., false), λ), vector(qb))
    return X, Q
end

function ∂g∂ʳposa(joint::Translational{T,N,N̄,Nl}, statea::State, stateb::State, λ, Δt) where {T,N,N̄,Nl}
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)
    ∂g∂ʳposa(joint, xa, qa, xb, qb, λ)
end

function ∂g∂ʳposa(joint::Translational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    X = -1.0 * transpose(rotation_matrix(qa))
    pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    ca_a = rotation_matrix(inv(qa)) * (xa) # body a com
    capb_a = pb_a - ca_a
    Q = - 1.0 * transpose(skew(capb_a))


    return [
            constraintmat(joint) * [X Q]
            -nullspacemat(joint) * [X Q]
            nullspacemat(joint) * [X Q]
           ]
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,Nl}, statea::State, stateb::State, λ, Δt) where {T,N,N̄,Nl}
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)
    ∂g∂ʳposb(joint, xa, qa, xb, qb, λ)
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    X = transpose(rotation_matrix(qa))
    pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    cbpb_b = pb_b - cb_b
    Q = transpose(skew(cbpb_b) * rotation_matrix(inv(qb) * qa))
    return [
            constraintmat(joint) * [X Q]
            -nullspacemat(joint) * [X Q]
            nullspacemat(joint) * [X Q]
           ]
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,Nl}, stateb::State, λ, Δt) where {T,N,N̄,Nl}
    xb, qb = posargs2(stateb)
    ∂g∂ʳposb(joint, xb, qb, λ)
end

function ∂g∂ʳposb(joint::Translational{T,N,N̄,Nl}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    X = transpose(I(3))
    pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    cbpb_b = pb_b - cb_b
    Q = transpose(skew(cbpb_b) * rotation_matrix(inv(qb)))
    return [
            constraintmat(joint) * [X Q]
            -nullspacemat(joint) * [X Q]
            nullspacemat(joint) * [X Q]
           ]
end

@inline function ∂g∂ʳself(joint::Translational{T,N,N̄,0}, λ) where {T,N,N̄}
    return Diagonal(1e-10 * sones(T,N))
end

@inline function ∂g∂ʳself(joint::Translational{T,N,N̄,Nl}, λ) where {T,N,N̄,Nl}
    [
     zeros(N, N + 4Nl);
     zeros(2Nl, N) Diagonal(ones(2Nl)) zeros(2Nl, 2Nl);
    ]
end

### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::Translational{T}, statea::State, stateb::State, Δt::T, clear::Bool) where T
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)

    Faw, τaa, Fbw, τbb = applyFτ(joint, joint.Fτ, xa, qa, xb, qb)
    statea.F2[end] += Faw
    statea.τ2[end] += τaa/2
    stateb.F2[end] += Fbw
    stateb.τ2[end] += τbb/2
    clear && (joint.Fτ = szeros(T,3))
    return
end
@inline function applyFτ(joint::Translational{T}, F::AbstractVector, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    vertices = joint.vertices

    Faw = vrotate(-F, qa) # in the world frame
    Fbw = -Faw # in the world frame
    Faa = vrotate(Faw, inv(qa)) # in local frame
    Fbb = vrotate(Fbw, inv(qb)) # in local frame

    pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point in b frame
    ca_a = rotation_matrix(inv(qa)) * (xa) # body a com in a frame
    ra = pb_a - ca_a
    τaa = torqueFromForce(Faa, ra) # in local coordinates
    τbb = torqueFromForce(Fbb, vertices[2]) # in local coordinates
    return Faw, τaa, Fbw, τbb
end

@inline function applyFτ!(joint::Translational{T}, stateb::State, Δt::T, clear::Bool) where T
    xb, qb = posargs2(stateb)

    Fbw, τbb = applyFτ(joint, joint.Fτ, xb, qb)
    stateb.F2[end] += Fbw
    stateb.τ2[end] += τbb/2
    clear && (joint.Fτ = szeros(T,3))
    return
end

@inline function applyFτ(joint::Translational{T}, F::AbstractVector, xb::AbstractVector, qb::UnitQuaternion) where T
    vertices = joint.vertices

    Fbw = F # in world frame
    Fbb = vrotate(Fbw, inv(qb)) # in b frame

    τbb = torqueFromForce(Fbb, vertices[2]) # in local coordinates
    return Fbw, τbb
end


## Forcing derivatives (for linearization)
# Control derivatives
@inline function ∂Fτ∂ua(joint::Translational, statea::State, stateb::State, Δt::T) where T
    vertices = joint.vertices
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)


    BFa = FiniteDiff.finite_difference_jacobian(F -> applyFτ(joint, F, xa, qa, xb, qb)[1], joint.Fτ)
    Bτa = 0.5 * FiniteDiff.finite_difference_jacobian(F -> applyFτ(joint, F, xa, qa, xb, qb)[2], joint.Fτ)

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Translational, statea::State, stateb::State, Δt::T) where T
    vertices = joint.vertices
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)

    BFb = FiniteDiff.finite_difference_jacobian(F -> applyFτ(joint, F, xa, qa, xb, qb)[3], joint.Fτ)
    Bτb = 0.5 * FiniteDiff.finite_difference_jacobian(F -> applyFτ(joint, F, xa, qa, xb, qb)[4], joint.Fτ)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Translational, stateb::State, Δt::T) where T
    vertices = joint.vertices
    xb, qb = posargs2(stateb)

    BFb = FiniteDiff.finite_difference_jacobian(F -> applyFτ(joint, F, xb, qb)[1], joint.Fτ)
    Bτb = 0.5 * FiniteDiff.finite_difference_jacobian(F -> applyFτ(joint, F, xb, qb)[2], joint.Fτ)

    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::Translational{T}, statea::State, stateb::State, Δt::T) where T
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXa = FiniteDiff.finite_difference_jacobian(xa -> applyFτ(joint, F, xa, qa, xb, qb)[1], xa)
    FaQa = FiniteDiff.finite_difference_jacobian(qa -> applyFτ(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[1], [qa.w, qa.x, qa.y, qa.z])
    τaXa = 0.5 * FiniteDiff.finite_difference_jacobian(xa -> applyFτ(joint, F, xa, qa, xb, qb)[2], xa)
    τaQa = 0.5 * FiniteDiff.finite_difference_jacobian(qa -> applyFτ(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[2], [qa.w, qa.x, qa.y, qa.z])
    FbXa = FiniteDiff.finite_difference_jacobian(xa -> applyFτ(joint, F, xa, qa, xb, qb)[3], xa)
    FbQa = FiniteDiff.finite_difference_jacobian(qa -> applyFτ(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[3], [qa.w, qa.x, qa.y, qa.z])
    τbXa = 0.5 * FiniteDiff.finite_difference_jacobian(xa -> applyFτ(joint, F, xa, qa, xb, qb)[4], xa)
    τbQa = 0.5 * FiniteDiff.finite_difference_jacobian(qa -> applyFτ(joint, F, xa, UnitQuaternion(qa..., false), xb, qb)[4], [qa.w, qa.x, qa.y, qa.z])

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::Translational{T}, statea::State, stateb::State, Δt::T) where T
    xa, qa = posargs2(statea)
    xb, qb = posargs2(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = FiniteDiff.finite_difference_jacobian(xb -> applyFτ(joint, F, xa, qa, xb, qb)[1], xb)
    FaQb = FiniteDiff.finite_difference_jacobian(qb -> applyFτ(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[1], [qb.w, qb.x, qb.y, qb.z])
    τaXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> applyFτ(joint, F, xa, qa, xb, qb)[2], xb)
    τaQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> applyFτ(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[2], [qb.w, qb.x, qb.y, qb.z])
    FbXb = FiniteDiff.finite_difference_jacobian(xb -> applyFτ(joint, F, xa, qa, xb, qb)[3], xb)
    FbQb = FiniteDiff.finite_difference_jacobian(qb -> applyFτ(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[3], [qb.w, qb.x, qb.y, qb.z])
    τbXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> applyFτ(joint, F, xa, qa, xb, qb)[4], xb)
    τbQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> applyFτ(joint, F, xa, qa, xb, UnitQuaternion(qb..., false))[4], [qb.w, qb.x, qb.y, qb.z])

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::Translational{T}, stateb::State, Δt::T) where T
    xb, qb = posargs2(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,4)
    τaXb = 0.5 * szeros(T,3,3)
    τaQb = 0.5 * szeros(T,3,4)
    FbXb = FiniteDiff.finite_difference_jacobian(xb -> applyFτ(joint, F, xb, qb)[1], xb)
    FbQb = FiniteDiff.finite_difference_jacobian(qb -> applyFτ(joint, F, xb, UnitQuaternion(qb..., false))[1], [qb.w, qb.x, qb.y, qb.z])
    τbXb = 0.5 * FiniteDiff.finite_difference_jacobian(xb -> applyFτ(joint, F, xb, qb)[2], xb)
    τbQb = 0.5 * FiniteDiff.finite_difference_jacobian(qb -> applyFτ(joint, F, xb, UnitQuaternion(qb..., false))[2], [qb.w, qb.x, qb.y, qb.z])

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end

## Position and velocity offsets
@inline function getPositionDelta(joint::Translational, body1::Component, body2::Component, x::SVector)
    Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational, body1::Component, body2::Component, v::SVector)
    Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
    return Δv
end

## Minimal coordinate calculation
@inline function minimalCoordinates(joint::Translational, body1::Component, body2::Component)
    statea = body1.state
    stateb = body2.state
    return minimalCoordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end
@inline function minimalCoordinates(joint::Translational, body1::Origin, body2::Component)
    stateb = body2.state
    return minimalCoordinates(joint, stateb.x2[1], stateb.q2[1])
end
@inline function minimalCoordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspacemat(joint) * _g(joint, xa, qa, xb, qb)
end
@inline function minimalCoordinates(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    return nullspacemat(joint) * _g(joint, xb, qb)
end
@inline function minimalVelocities(joint::Translational, body1::Component, body2::Component)
    statea = body1.state
    stateb = body2.state
    return minimalVelocities(joint, statea.x2[1], statea.q2[1], statea.v15, statea.ϕ15, stateb.x2[1], stateb.q2[1], stateb.v15, stateb.ϕ15)
end

@inline function minimalVelocities(joint::Translational, body1::Origin, body2::Component)
    stateb = body2.state
    return minimalVelocities(joint, stateb.q2[1], stateb.v15, stateb.ϕ15)
end

@inline function minimalVelocities(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector)
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa))
    Δv = vrotate(Δvw, inv(qa))
    return nullspacemat(joint) * Δv
end

@inline function minimalVelocities(joint::Translational, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector)
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb)
    Δv = Δvw
    return nullspacemat(joint) * Δv
end
