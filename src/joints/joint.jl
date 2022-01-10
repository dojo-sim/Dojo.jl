abstract type Joint{T,Nλ,Nb,N} end

### General functions
getT(joint::Joint{T}) where T = T
Base.length(joint::Joint{T,Nλ}) where {T,Nλ} = Nλ
Base.zero(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, Nλ, 6)
@inline g(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, Nλ)

λlength(joint::Joint{T,Nλ}) where {T,Nλ} = Nλ
blength(joint::Joint{T,Nλ,Nb}) where {T,Nλ,Nb} = Nb
ηlength(joint::Joint{T,Nλ,Nb,N}) where {T,Nλ,Nb,N} = N

function get_sγ(joint::Joint{T,Nλ,Nb}, η) where {T,Nλ,Nb}
    s = η[SVector{Nb,Int}(1:Nb)]
    γ = η[SVector{Nb,Int}(Nb .+ (1:Nb))]
    return s, γ
end

function λindex(joint::Joint{T,Nλ,Nb,N}, s::Int) where {T,Nλ,Nb,N}
    ind = SVector{N,Int}(s+1:s+N)
    return ind
end

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳposa(joint::Joint, body1::Body, body2::Body, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂ʳposa(joint, body1.state, body2.state, λ, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Joint, body1::Body, body2::Body, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂ʳposb(joint, body1.state, body2.state, λ, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳposb(joint::Joint, body1::Origin, body2::Body, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂ʳposb(joint, body2.state, λ, Δt)
    else
        return zero(joint)
    end
end

## Discrete-time velocity derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳvela(joint::Joint, body1::Body, body2::Body, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂ʳvela(joint, body1.state, body2.state, λ, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Joint, body1::Body, body2::Body, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂ʳvelb(joint, body1.state, body2.state, λ, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::Joint, body1::Origin, body2::Body, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂ʳvelb(joint, body2.state, λ, Δt)
    else
        return zero(joint)
    end
end

# @inline function offdiagonal∂damper∂ʳvel(joint::Joint, body1::Body, body2::Body, childid)
#     if body2.id == childid
#         return offdiagonal∂damper∂ʳvel(joint, body1.state, body2.state)
#     else
#         return zero(body2)
#     end
# end
# @inline function offdiagonal∂damper∂ʳvel(joint::Joint, body1::Origin, body2::Body, childid)
#     if body2.id == childid
#         return offdiagonal∂damper∂ʳvel(joint, body2.state)
#     else
#         return zero(body2)
#     end
# end

### Springs and Dampers (for dynamics)
## Wrappers 1
@inline function springforcea(joint::Joint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return springforcea(joint, body1, body2, Δt)
    else
        return springforce(joint)
    end
end
@inline function springforceb(joint::Joint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return springforceb(joint, body1, body2, Δt)
    else
        return springforce(joint)
    end
end
@inline function springforceb(joint::Joint, body1::Origin, body2::Body, Δt, childid)
    if body2.id == childid
        return springforceb(joint, body2, Δt)
    else
        return springforce(joint)
    end
end

@inline function damperforcea(joint::Joint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return damperforcea(joint, body1, body2, Δt)
    else
        return damperforce(joint)
    end
end
@inline function damperforceb(joint::Joint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return damperforceb(joint, body1, body2, Δt)
    else
        return damperforce(joint)
    end
end
@inline function damperforceb(joint::Joint, body1::Origin, body2::Body, Δt, childid)
    if body2.id == childid
        return damperforceb(joint, body2, Δt)
    else
        return damperforce(joint)
    end
end

## Wrappers 2
@inline springforce(joint::Joint{T}) where {T} = szeros(T, 6) # TODO zero function?

@inline damperforce(joint::Joint{T}) where {T} = szeros(T, 6) # TODO zero function?

# ### Forcing (for dynamics)
## Wrappers
@inline function applyFτ!(joint::Joint, body1::Body, body2::Body, Δt::T, clear::Bool) where T
    applyFτ!(joint, body1.state, body2.state, Δt, clear)
    return
end

@inline function applyFτ!(joint::Joint, ::Origin, body2::Body, Δt::T, clear::Bool) where T
    applyFτ!(joint, body2.state, Δt, clear)
    return
end

Joint0 = Joint{T,0} where T
Joint1 = Joint{T,1} where T
Joint2 = Joint{T,2} where T
Joint3 = Joint{T,3} where T

# Base.show(io::IO, joint::Joint) = summary(io, joint)

### Constaint and nullspace matrices
@inline constraintmat(::Joint0{T}) where T = szeros(T,0,3)
@inline nullspacemat(::Joint0{T}) where T = SMatrix{3,3,T,9}(I)
@inline constraintmat(joint::Joint1) = joint.V3
@inline nullspacemat(joint::Joint1) = joint.V12
@inline constraintmat(joint::Joint2) = joint.V12
@inline nullspacemat(joint::Joint2) = joint.V3
@inline constraintmat(::Joint3{T}) where T = SMatrix{3,3,T,9}(I)
@inline nullspacemat(::Joint3{T}) where T = szeros(T,0,3)

### Constraints and derivatives
## Position level constraint wrappers
@inline g(joint::Joint, body1::Body, body2::Body, λ, Δt) = g(joint, body1.state, body2.state, λ, Δt)
@inline g(joint::Joint, body1::Origin, body2::Body, λ, Δt) = g(joint, body2.state, λ, Δt)

### Constraints and derivatives
## Discrete-time position wrappers (for dynamics)
g(joint::Joint, statea::State, stateb::State, λ, Δt) = g(joint, posargs3(statea, Δt)..., posargs3(stateb, Δt)..., λ)
g(joint::Joint, stateb::State, λ, Δt) = g(joint, posargs3(stateb, Δt)..., λ)

@inline function ∂g∂ʳself(joint::Joint{T,Nλ}, λ) where {T,Nλ}
    return Diagonal(+1.00e-10 * sones(T,Nλ))
end

## Discrete-time position derivatives (for dynamics)
# Wrappers 1

@inline ∂g∂posa(joint::Joint, body1::Body, body2::Body, λ, Δt) = ∂g∂posa(joint, posargs3(body1.state, Δt)..., posargs3(body2.state, Δt)..., λ)
@inline ∂g∂posb(joint::Joint, body1::Body, body2::Body, λ, Δt) = ∂g∂posb(joint, posargs3(body1.state, Δt)..., posargs3(body2.state, Δt)..., λ)
@inline ∂g∂posb(joint::Joint, body1::Origin, body2::Body, λ, Δt) = ∂g∂posb(joint, posargs3(body2.state, Δt)..., λ)

# Wrappers 2
∂g∂ʳposa(joint::Joint, statea::State, stateb::State, λ, Δt) = ∂g∂ʳposa(joint, posargs2(statea)..., posargs2(stateb)..., λ)
∂g∂ʳposb(joint::Joint, statea::State, stateb::State, λ, Δt) = ∂g∂ʳposb(joint, posargs2(statea)..., posargs2(stateb)..., λ)
∂g∂ʳposb(joint::Joint, stateb::State, λ, Δt) = ∂g∂ʳposb(joint, posargs2(stateb)..., λ)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳposa(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ)
    X, Q = ∂g∂posa(joint, xa, qa, xb, qb, λ)
    Q = Q * LVᵀmat(qa)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Joint, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ)
    X, Q = ∂g∂posb(joint, xa, qa, xb, qb, λ)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Joint, xb::AbstractVector, qb::UnitQuaternion, λ)
    X, Q = ∂g∂posb(joint, xb, qb, λ)
    Q = Q * LVᵀmat(qb)

    return [X Q]
end

# Wrappers 2
∂g∂ʳvela(joint::Joint, statea::State, stateb::State, λ, Δt) = ∂g∂ʳvela(joint, posargs3(statea, Δt)..., posargs3(stateb, Δt)..., fullargssol(statea)..., λ, Δt)
∂g∂ʳvelb(joint::Joint, statea::State, stateb::State, λ, Δt) = ∂g∂ʳvelb(joint, posargs3(statea, Δt)..., posargs3(stateb, Δt)..., fullargssol(stateb)..., λ, Δt)
∂g∂ʳvelb(joint::Joint, stateb::State, λ, Δt) = ∂g∂ʳvelb(joint, posargs3(stateb, Δt)..., fullargssol(stateb)..., λ, Δt)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvela(joint::Joint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1a::AbstractVector, v1a::AbstractVector, q1a::UnitQuaternion, ω1a::AbstractVector, λ, Δt,
    )

    X, Q = ∂g∂posa(joint, x2a, q2a, x2b, q2b, λ)
    V = X * Δt
    # Ω = Q * Lmat(q1a) * derivωbar(ω1a, Δt) * Δt / 2
    Ω = Q * Lmat(q1a) * derivcayley(ω1a)

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Joint, x2a::AbstractVector, q2a::UnitQuaternion, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, λ, Δt
    )

    X, Q = ∂g∂posb(joint, x2a, q2a, x2b, q2b, λ)
    V = X * Δt
    # Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt) * Δt / 2
    Ω = Q * Lmat(q1b) * derivcayley(ω1b)

    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Joint, x2b::AbstractVector, q2b::UnitQuaternion,
        x1b::AbstractVector, v1b::AbstractVector, q1b::UnitQuaternion, ω1b::AbstractVector, λ, Δt
    )

    X, Q = ∂g∂posb(joint, x2b, q2b, λ)
    V = X * Δt
    # Ω = Q * Lmat(q1b) * derivωbar(ω1b, Δt) * Δt / 2
    Ω = Q * Lmat(q1b) * derivcayley(ω1b)

    return [V Ω]
end

### Force derivatives (for linearization)
## Forcing
@inline function setForce!(joint::Joint, Fτ::SVector)
    joint.Fτ = zerodimstaticadjoint(nullspacemat(joint)) * Fτ
    return
end
@inline setForce!(joint::Joint) = return


@inline function addForce!(joint::Joint, Fτ::SVector)
    joint.Fτ += zerodimstaticadjoint(nullspacemat(joint)) * Fτ
    return
end
@inline addForce!(joint::Joint) = return

## Derivative wrappers
@inline function ∂Fτ∂ua(joint::Joint, body1::Body, body2::Body, Δt, childid)
    return ∂Fτ∂ua(joint, body1.state, body2.state, Δt) * zerodimstaticadjoint(nullspacemat(joint))
end
@inline function ∂Fτ∂ub(joint::Joint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return ∂Fτ∂ub(joint, body1.state, body2.state, Δt) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Joint, body1::Origin, body2::Body, Δt, childid)
    if body2.id == childid
        return return ∂Fτ∂ub(joint, body2.state, Δt) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline ∂Fτ∂ub(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, 6, 3 - Nλ) # TODO zero function?


### Minimal coordinates
@inline minimalCoordinates(joint::Joint{T,Nλ}) where {T,Nλ} = szeros(T, 3 - Nλ)

## Limits
function add_limits(mech::Mechanism, eq::EqualityConstraint;
    # NOTE: this only works for joints between serial chains (ie, single child joints)
    tra_limits=eq.constraints[1].joint_limits,
    rot_limits=eq.constraints[1].joint_limits)

    # update translational
    tra = eq.constraints[1]
    T = typeof(tra).parameters[1]
    Nλ = typeof(tra).parameters[2]
    Nb½ = length(tra_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    tra_limit = (Translational{T,Nλ,Nb,N,Nb½,N̄λ}(tra.V3, tra.V12, tra.vertices, tra.spring, tra.damper, tra.spring_offset, tra_limits, tra.spring_type, tra.Fτ), eq.parentid, eq.childids[1])

    # update rotational
    rot = eq.constraints[2]
    T = typeof(rot).parameters[1]
    Nλ = typeof(rot).parameters[2]
    Nb½ = length(rot_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    rot_limit = (Rotational{T,Nλ,Nb,N,Nb½,N̄λ}(rot.V3, rot.V12, rot.qoffset, rot.spring, rot.damper, rot.spring_offset, rot_limits, rot.spring_type, rot.Fτ), eq.parentid, eq.childids[1])
    EqualityConstraint((tra_limit, rot_limit); name=eq.name)
end
