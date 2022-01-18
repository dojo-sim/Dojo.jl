abstract type Joint{T,Nλ,Nb,N} end

## General functions
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
@inline function Ga(joint::Joint, body1::Component, body2::Component, childid, λ, Δt)
    if body2.id == childid
        return Ga(joint, posargs2(body1.state)..., posargs2(body2.state)..., λ)
    else
        return zero(joint)
    end
end

@inline function Gb(joint::Joint, body1::Component, body2::Component, childid, λ, Δt)
    if body2.id == childid
        return Gb(joint, posargs2(body1.state)..., posargs2(body2.state)..., λ)
    else
        return zero(joint)
    end
end

## Discrete-time velocity derivatives (for dynamics)
@inline function ∂g∂a(joint::Joint, body1::Component, body2::Component, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂a(joint, posargs3(body1.state, Δt)..., posargs3(body2.state, Δt)..., λ)
    else
        return zero(joint)
    end
end
@inline function ∂g∂b(joint::Joint, body1::Component, body2::Component, childid, λ, Δt)
    if body2.id == childid
        return ∂g∂b(joint, posargs3(body1.state, Δt)..., posargs3(body2.state, Δt)..., λ)
        
    else
        return zero(joint)
    end
end


### Springs and Dampers (for dynamics)
@inline function springforcea(joint::Joint, body1::Component, body2::Component, Δt, childid)
    if body2.id == childid
        return springforcea(joint, body1, body2, Δt)
    else
        return szeros(T, 6)
    end
end
@inline function springforceb(joint::Joint, body1::Component, body2::Component, Δt, childid)
    if body2.id == childid
        return springforceb(joint, body1, body2, Δt)
    else
        return szeros(T, 6)
    end
end

@inline function damperforcea(joint::Joint, body1::Component, body2::Component, Δt, childid)
    if body2.id == childid
        return damperforcea(joint, body1, body2, Δt)
    else
        return szeros(T, 6)
    end
end

@inline function damperforceb(joint::Joint, body1::Component, body2::Component, Δt, childid)
    if body2.id == childid
        return damperforceb(joint, body1, body2, Δt)
    else
        return szeros(T, 6)
    end
end

# ### Forcing (for dynamics)
@inline function applyFτ!(joint::Joint, body1::Component, body2::Component, Δt::T, clear::Bool) where T
    applyFτ!(joint, body1.state, body2.state, Δt, clear)
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
@inline g(joint::Joint, body1::Component, body2::Component, λ, Δt) = g(joint, body1.state, body2.state, λ, Δt)

### Constraints and derivatives
## Discrete-time position wrappers (for dynamics)
g(joint::Joint, statea::State, stateb::State, λ, Δt) = g(joint, posargs3(statea, Δt)..., posargs3(stateb, Δt)..., λ)

@inline function ∂g∂z(joint::Joint{T,Nλ}, λ) where {T,Nλ}
    return Diagonal(+1.00e-10 * sones(T,Nλ))
end

## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline ∂g∂a(joint::Joint, body1::Component, body2::Component, λ, Δt) = ∂g∂a(joint, posargs3(body1.state, Δt)..., posargs3(body2.state, Δt)..., λ)
@inline ∂g∂b(joint::Joint, body1::Component, body2::Component, λ, Δt) = ∂g∂b(joint, posargs3(body1.state, Δt)..., posargs3(body2.state, Δt)..., λ)

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
@inline function ∂Fτ∂ua(joint::Joint, body1::Component, body2::Component, Δt, childid)
    return ∂Fτ∂ua(joint, body1.state, body2.state, Δt) * zerodimstaticadjoint(nullspacemat(joint))
end

@inline function ∂Fτ∂ub(joint::Joint{T,Nλ}, body1::Component, body2::Component, Δt, childid) where {T,Nλ}
    if body2.id == childid
        return ∂Fτ∂ub(joint, body1.state, body2.state, Δt) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return szeros(T, 6, 3 - Nλ)
    end
end

## Minimal coordinates
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
