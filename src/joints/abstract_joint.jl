abstract type AbstractJoint{T,N} end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::AbstractJoint{T,N}) where {T,N}
    summary(io, constraint)
end

### General functions
getT(joint::AbstractJoint{T}) where T = T
Base.length(joint::AbstractJoint{T,N}) where {T,N} = N
Base.zero(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N, 6)
@inline g(joint::AbstractJoint{T,N}) where {T,N} = szeros(T, N)


## Discrete-time position derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳposa(joint::AbstractJoint, body1::Body, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳposa(joint, body1.state, body2.state, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, body1::Body, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳposb(joint, body1.state, body2.state, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳposb(joint::AbstractJoint, body1::Origin, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳposb(joint, body2.state, Δt)
    else
        return zero(joint)
    end
end




## Discrete-time velocity derivatives (for dynamics)
# Wrappers 1
@inline function ∂g∂ʳvela(joint::AbstractJoint, body1::Body, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳvela(joint, body1.state, body2.state, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, body1::Body, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳvelb(joint, body1.state, body2.state, Δt)
    else
        return zero(joint)
    end
end
@inline function ∂g∂ʳvelb(joint::AbstractJoint, body1::Origin, body2::Body, childid, Δt)
    if body2.id == childid
        return constraintmat(joint) * ∂g∂ʳvelb(joint, body2.state, Δt)
    else
        return zero(joint)
    end
end







@inline function offdiagonal∂damper∂ʳvel(joint::AbstractJoint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return offdiagonal∂damper∂ʳvel(joint, body1.state, body2.state)
    else
        return zero(body2)
    end
end
@inline function offdiagonal∂damper∂ʳvel(joint::AbstractJoint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return offdiagonal∂damper∂ʳvel(joint, body2.state)
    else
        return zero(body2)
    end
end



### Springs and Dampers (for dynamics)
## Wrappers 1
@inline function springforcea(joint::AbstractJoint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return springforcea(joint, body1.state, body2.state, Δt)
    else
        return springforce(joint)
    end
end
@inline function springforceb(joint::AbstractJoint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return springforceb(joint, body1.state, body2.state, Δt)
    else
        return springforce(joint)
    end
end
@inline function springforceb(joint::AbstractJoint, body1::Origin, body2::Body, Δt, childid)
    if body2.id == childid
        return springforceb(joint, body2.state, Δt)
    else
        return springforce(joint)
    end
end

@inline function damperforcea(joint::AbstractJoint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return damperforcea(joint, body1.state, body2.state, Δt)
    else
        return damperforce(joint)
    end
end
@inline function damperforceb(joint::AbstractJoint, body1::Body, body2::Body, Δt, childid)
    if body2.id == childid
        return damperforceb(joint, body1.state, body2.state, Δt)
    else
        return damperforce(joint)
    end
end
@inline function damperforceb(joint::AbstractJoint, body1::Origin, body2::Body, Δt, childid)
    if body2.id == childid
        return damperforceb(joint, body2.state, Δt)
    else
        return damperforce(joint)
    end
end

## Wrappers 2
@inline springforce(joint::AbstractJoint{T}) where {T} = szeros(T, 6) # TODO zero function?
# springforcea(joint::AbstractJoint, statea::State, stateb::State) = springforcea(joint, posargsk(statea)..., posargsk(stateb)...)
# springforceb(joint::AbstractJoint, statea::State, stateb::State) = springforceb(joint, posargsk(statea)..., posargsk(stateb)...)
# springforceb(joint::AbstractJoint, stateb::State) = springforceb(joint, posargsk(stateb)...)

@inline damperforce(joint::AbstractJoint{T}) where {T} = szeros(T, 6) # TODO zero function?
# damperforcea(joint::AbstractJoint, statea::State, stateb::State) = damperforcea(joint, fullargssol(statea)..., fullargssol(stateb)...)
# damperforceb(joint::AbstractJoint, statea::State, stateb::State) = damperforceb(joint, fullargssol(statea)..., fullargssol(stateb)...)
# damperforceb(joint::AbstractJoint, stateb::State) = damperforceb(joint, fullargssol(stateb)...)


### Forcing (for dynamics)
## Wrappers
@inline function applyFτ!(joint::AbstractJoint, body1::Body, body2::Body, Δt::T, clear::Bool) where T
    applyFτ!(joint, body1.state, body2.state, Δt, clear)
    return
end
@inline function applyFτ!(joint::AbstractJoint, ::Origin, body2::Body, Δt::T, clear::Bool) where T
    applyFτ!(joint, body2.state, Δt, clear)
    return
end
