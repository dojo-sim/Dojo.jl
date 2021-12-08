
# ################################################################################
# # Reward
# ################################################################################
#
# """
#      Reward constructor. Provides a simple way to construct a reward function
#      for each envirionment: atlas, snake, dice, etc.
# """
# function getreward(model::Symbol; kwargs...)
#     reward_fct = eval(Symbol(:getreward, model))(; kwargs...)
#     return reward_fct
# end
#
# function getrewarddice(;)
#     reward_fct(s, a) = 0.0
#     return reward_fct
# end

################################################################################
# Environment
################################################################################

abstract type Environment{T,M,A,O} end

function make(model::String; kwargs...)
    return eval(Symbol(model))(; kwargs...)
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
# Pendulum
################################################################################

mutable struct Pendulum{T,M,A,O} <: Environment{T,M,A,O}
    mechanism::M
    aspace::A
    ospace::O
    x::AbstractVector{T}
    last_u::AbstractVector{T}
    nx::Int
    nu::Int
    no::Int
    max_speed::T # never used because we do not clip the pendulum angular velocity
    max_torque::T
    dt::T
    g::T
    m::T
    l::T
    rng::MersenneTwister
    vis::Visualizer
end

function Pendulum(; max_speed::T = 8.0, max_torque::T = 40.0,
        dt::T = 0.05, g::T = -10.0, m::T = 1.0, l::T = 1.0, s::Int = 1, vis::Visualizer = Visualizer()) where {T}
    mechanism = getmechanism(:pendulum, Δt = dt, g = g, m = m, l = l, damper = 5.0)
    nx = minCoordDim(mechanism)
    nu = controldim(mechanism)
    no = 3

    high = [1.0, 1.0, max_speed]
    aspace = BoxSpace(controldim(mechanism), low = [-max_torque], high = [max_torque])
    ospace = BoxSpace(no, low = -high, high = high)
    rng = MersenneTwister(s)
    x = Inf * ones(nx)
    last_u = Inf * ones(nu)
    build_robot(vis, mechanism)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = Pendulum{TYPES...}(mechanism, aspace, ospace, x, last_u, nx, nu, no,
        max_speed, max_torque, dt, g, m, l, rng, vis)
    seed(env, s = s)
    return env
end

function seed(env::Pendulum{T}; s = 0) where {T}
    env.rng = MersenneTwister(s)
    return nothing
end

function reset(env::Pendulum{T}; x = nothing) where {T}
    if x != nothing
        env.x = x
    else
        high = [π, 1.0]
        low = -high
        env.x = rand(env.rng, env.nx) .* (high .- low) .+ low
        env.last_u .= Inf
    end
    return _get_obs(env)
end

function _get_obs(env::Pendulum{T}) where {T}
    θ, ω = env.x
    return [cos(θ), sin(θ), ω]
end

function step(env::Pendulum{T}, u::AbstractVector{T}) where {T}
    mechanism = env.mechanism
    Δt = mechanism.Δt
    x0 = env.x


    u0 = clamp.(u, -env.max_torque, env.max_torque)
    env.last_u = u0  # for rendering

    z0 = min2max(mechanism, x0)
    z1 = simon_step!(mechanism, z0, Δt * u0; ϵ = 1e-6, newtonIter = 100,
        lineIter = 10, verbose = false, btol = 1e-6, undercut = Inf)
    env.x = max2min(mechanism, z1)

    # Compute cost function
    θ0, ω0 = x0
    costs = angle_normalize(θ0 - π)^2 + 1e-1 * ω0^2 + 1e-3 * u[1]^2 # angle_normalize enforces angle ∈ [-π, π]

    info = Dict()
    return _get_obs(env), -costs, false, info
end

function angle_normalize(x)
    return ((x + π) % (2 * π)) - π
end

function render(env::Pendulum, mode="human")
    z = min2max(env.mechanism, env.x)
    set_robot(env.vis, env.mechanism, z)
    return nothing
end

function close(env::Pendulum{M}; kwargs...) where {M}
    return nothing
end



# ################################################################################
# # Environment
# ################################################################################
#
# mutable struct Environment{M,A,O}
#     mechanism::M
#     model::Symbol
#     reward_fct::Any # reward function R(s, a) = r
#     iter::Int # current time step count
#     max_iter::Int # number of time steps per episode
#     aspace::A
#     ospace::O
#     # vis::Visualizer
# end
#
#
# function Environment(mechanism::Mechanism, model::Symbol, reward_fct::Any; max_iter::Int = 100)
#     iter = 0
#     aspace = BoxSpace(controldim(mechanism))
#     ospace = BoxSpace(maxCoordDim(mechanism))
#     TYPE = [typeof(mechanism), typeof(aspace), typeof(ospace)]
#     env = Environment{TYPE...}(mechanism, model, reward_fct, iter, max_iter, aspace, ospace)
#     return env
# end
#
# function make(model::String; kwargs...)
#     model = Symbol(model)
#     mechanism = getmechanism(model; kwargs...)
#     reward_fct = getreward(model)
#     env = Environment(mechanism, model, reward_fct)
#     reset(env)
#     return env
# end
#
# function reset(env::Environment{M}; kwargs...) where {M}
#     env.iter = 0
#     initialize!(env.mechanism, Symbol(env.model); kwargs...)
#     observation = getMaxState(env.mechanism)
#     return observation
# end
#
# function render(env::Environment{M}; kwargs...) where {M}
#
#     return nothing
# end
#
# function close(env::Environment{M}; kwargs...) where {M}
#     # visualizer stuff
#     return nothing
# end
#
# function need_to_reset(env::Environment)
#     done = env.iter >= env.max_iter - 1
#     return done
# end
#
# function step(env::Environment{M}, action) where {M}
#     env.iter += 1
#     mechanism = env.mechanism
#     state = getMaxState(mechanism)
#     observation = getMaxState(mechanism)
#
#     reward = env.reward_fct(state, action)
#     done = need_to_reset(env)
#     info = Dict()
#     return observation, reward, done, info
# end
