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

struct BoxSpace{T,N} <: Space{T,N}
    n::Int # box dimension
    low::Vector{T} # minimum value
    high::Vector{T} # maximum value
end

function BoxSpace(n::Int; low::Vector{T}=-ones(n), high::Vector{T}=ones(n)) where T
    return BoxSpace{T,n}(n, low, high)
end

function sample(s::BoxSpace{T,N}) where {T,N}
    return rand(T, N) .* (s.high .- s.low) .+ s.low
end

function contains(s::BoxSpace{T,N}, v::Vector{T}) where {T,N}
    all(v .>= s.low) && all(v .<= s.high)
end

################################################################################
# Step
################################################################################
step(env::Environment, u) = step(env, env.x, u)

function f(y, env::Environment, x, u, w)
	y .= step(env, x, u)[1]
end

function fx(fx, env::Environment, x, u, w)
	step(env, x, u, diff=true)
    fx .= env.fx
end

function fu(fu, env::Environment, x, u, w)
	step(env, x, u, diff=true)
	fu .= env.fu
end

################################################################################
# Environments
################################################################################
include("pendulum/methods/env.jl")
include("hopper/methods/env.jl")













