
################################################################################
# Reward
################################################################################

"""
     Reward constructor. Provides a simple way to construct a reward function
     for each envirionment: atlas, snake, dice, etc.
"""
function getreward(model::Symbol; kwargs...)
    reward_fct = eval(Symbol(:getreward, model))(; kwargs...)
    return reward_fct
end

function getrewarddice(;)
    reward_fct(s, a) = 0.0
    return reward_fct
end


################################################################################
# Imports
################################################################################

import Base.contains
import Base.reset
import Base.step
import Dojo.MeshCat.render

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
# Environment
################################################################################

mutable struct Environment{M,A,O}
    mechanism::M
    model::Symbol
    reward_fct::Any # reward function R(s, a) = r
    iter::Int # current time step count
    max_iter::Int # number of time steps per episode
    aspace::A
    ospace::O
    # vis::Visualizer
end


function Environment(mechanism::Mechanism, model::Symbol, reward_fct::Any; max_iter::Int = 100)
    iter = 0
    aspace = BoxSpace(controldim(mechanism))
    ospace = BoxSpace(maxCoordDim(mechanism))
    TYPE = [typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = Environment{TYPE...}(mechanism, model, reward_fct, iter, max_iter, aspace, ospace)
    return env
end

function make(model::String; kwargs...)
    model = Symbol(model)
    mechanism = getmechanism(model; kwargs...)
    reward_fct = getreward(model)
    env = Environment(mechanism, model, reward_fct)
    reset(env)
    return env
end

function reset(env::Environment{M}; kwargs...) where {M}
    env.iter = 0
    initialize!(env.mechanism, Symbol(env.model); kwargs...)
    observation = getMaxState(env.mechanism)
    return observation
end

function render(env::Environment{M}; kwargs...) where {M}

    return nothing
end

function close(env::Environment{M}; kwargs...) where {M}
    # visualizer stuff
    return nothing
end

function need_to_reset(env::Environment)
    done = env.iter >= env.max_iter - 1
    return done
end

function step(env::Environment{M}, action) where {M}
    env.iter += 1
    mechanism = env.mechanism
    state = getMaxState(mechanism)
    observation = getMaxState(mechanism)

    reward = env.reward_fct(state, action)
    done = need_to_reset(env)
    info = Dict()
    return observation, reward, done, info
end


