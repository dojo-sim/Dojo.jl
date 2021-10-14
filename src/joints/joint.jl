# Standard translational or rotational joint
abstract type Joint{T,N} <: AbstractJoint{T,N} end

Joint0 = Joint{T,0} where T
Joint1 = Joint{T,1} where T
Joint2 = Joint{T,2} where T
Joint3 = Joint{T,3} where T

Base.show(io::IO, joint::Joint) = summary(io, joint)

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
## Discrete-time position wrappers (for dynamics)
g(joint::Joint, statea::State, stateb::State, Δt) = g(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
g(joint::Joint, stateb::State, Δt) = g(joint, posargsnext(stateb, Δt)...)
g(joint::Joint, statea::State, stateb::State) = g(joint, posargsc(statea)..., posargsc(stateb)...)
g(joint::Joint, stateb::State) = g(joint, posargsc(stateb)...)

@inline g(joint::Joint{T,N}) where {T,N} = szeros(T, N)


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
@inline function ∂Fτ∂ua(joint::Joint, body1::Body, body2::Body, childid)
    return ∂Fτ∂ua(joint, body1.state, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
end
@inline function ∂Fτ∂ub(joint::Joint, body1::Body, body2::Body, childid)
    if body2.id == childid
        return ∂Fτ∂ub(joint, body1.state, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end
end
@inline function ∂Fτ∂ub(joint::Joint, body1::Origin, body2::Body, childid)
    if body2.id == childid
        return return ∂Fτ∂ub(joint, body2.state) * zerodimstaticadjoint(nullspacemat(joint))
    else
        return ∂Fτ∂ub(joint)
    end
end

@inline ∂Fτ∂ub(joint::Joint{T,N}) where {T,N} = szeros(T, 6, 3 - N) # TODO zero function?


### Minimal coordinates
@inline minimalCoordinates(joint::Joint{T,N}) where {T,N} = szeros(T, 3 - N)







