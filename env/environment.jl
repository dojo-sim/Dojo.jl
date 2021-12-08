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

function make(model::String; kwargs...)
    return eval(Symbol(model))(; kwargs...)
end

struct Environment{X,T,M,A,O,I}
    mechanism::M
    mode::Symbol
    aspace::A
    ospace::O
    x::Vector{T}
    fx::Matrix{T}
    fu::Matrix{T}
    u_prev::Vector{T}
    nx::Int
    nu::Int
    no::Int
    info::I
    rng::Vector{MersenneTwister}
    vis::Visualizer
    opts_step::InteriorPointOptions{T}
    opts_grad::InteriorPointOptions{T}
end

function reset(env::Environment{X}; x=nothing) where X
    initialize!(X, env.mechanism)
    if x != nothing
        env.x = x
    else
        if env.mode == :min
            env.x .= getMinState(env.mechanism)
        elseif env.mode == :max
            env.x .= getMaxState(env.mechanism)
        end
        env.u_prev .= 0.0
    end
    return _get_obs(env)
end

function _get_obs(env::Environment)
    return env.x
end

function step(env::Environment, x, u; diff=false)
    mechanism = env.mechanism
    Δt = mechanism.Δt

    x0 = x
    env.u_prev .= u  # for rendering

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, env.control_mask' * u; opts=env.opts_step)
    env.x .= env.mode == :min ? max2min(mechanism, z1) : z1

    # Compute cost
    costs = reward(env, x, u)

    # Gradients
    if diff
        if env.mode == :min
            fx, fu = getMinGradients!(env.mechanism, z0, env.control_mask' * u, opts=env.opts_grad)
        elseif env.mode == :max
            fx, fu = getMaxGradients!(env.mechanism, z0, env.control_mask' * u, opts=env.opts_grad)
        end
        env.fx .= fx
        env.fu .= fu * env.control_mask'
    end

    info = Dict()

    return _get_obs(env), -costs, false, info
end

function render(env::Environment, mode="human")
    z = env.mode == :min ? min2max(env.mechanism, env.x) : env.x
    set_robot(env.vis, env.mechanism, z)
    return nothing
end

function seed(env::Environment; s=0)
    env.rng[1] = MersenneTwister(s)
    return nothing
end

reward(env::Environment, x, u) = 0.0

function close(env::Environment; kwargs...)
    return nothing
end

################################################################################
# Space
################################################################################

abstract type Space{T,N} end

mutable struct BoxSpace{T,N} <: Space{T,N,}
    n::Int # box dimension
    low::AbstractVector{T} # minimum value
    high::AbstractVector{T} # maximum value
    shape::Tuple{Int} # this is always (n,), it's needed to interface with Stable-Baselines
    dtype::DataType # this is always T, it's needed to interface with Stable-Baselines
end

function BoxSpace(n::Int; low::AbstractVector{T} = -ones(n), high::AbstractVector{T} = ones(n)) where {T}
    return BoxSpace{T,n}(n, low, high, (n,), T)
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

function f(y, env::Environment, x, u, w)
	y .= step(env, x, u)[1]
end

function fx(fx, env::Environment, x, u, w)
	step(env, x, u, diff=true)
    fx .= env.fx
end

function fu(fu, env::Environment, x, u, w)
	# step(env, x, u, diff=true) # this is run in fx
	fu .= env.fu
end

################################################################################
# Environments
# ##############################################################################
# include("pendulum/methods/env.jl")
# include("hopper/methods/env.jl")
# include("quadruped/methods/env.jl")
