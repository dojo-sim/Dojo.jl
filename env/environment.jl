################################################################################
# Imports
################################################################################

import Base.contains
import Base.reset
import Base.step
import Dojo.MeshCat.render

################################################################################
# Environment
################################################################################

abstract type Environment{T,M,A,O} end

function make(model::String; kwargs...)
    return eval(Symbol(model))(; kwargs...)
end

################################################################################
# Space
################################################################################

abstract type Space{T,N} end

mutable struct BoxSpace{T,N} <: Space{T,N}
    n::Int # box dimension
    low::AbstractVector{T} # minimum value
    high::AbstractVector{T} # maximum value
end

function BoxSpace(n::Int; low::AbstractVector{T} = -ones(n), high::AbstractVector{T} = ones(n)) where {T}
    return BoxSpace{T,n}(n, low, high)
end

function sample(s::BoxSpace{T,N}) where {T,N}
    return rand(T,N) .* (s.high .- s.low) .+ s.low
end

function contains(s::BoxSpace{T,N}, v::AbstractVector{T}) where {T,N}
    all(v .>= s.low) && all(v .<= s.high)
end

################################################################################
# Step
################################################################################
step(env::Environment, u) = step(env, env.x, u)

################################################################################
# Environments
################################################################################

include("pendulum/methods/env.jl")












