# include("environment.jl")


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

function getrewardpendulum(;)
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

abstract type Space43{T,N} end

mutable struct BoxSpace43{T,N} <: Space43{T,N}
    n::Int # box dimension
    low::AbstractVector{T} # minimum value
    high::AbstractVector{T} # maximum value
end

function BoxSpace43(n::Int; low::AbstractVector{T} = -ones(n), high::AbstractVector{T} = ones(n)) where {T}
    return BoxSpace43{T,n}(n, low, high)
end

function sample(s::BoxSpace43{T,N}) where {T,N}
    return rand(T,N) .* (s.high .- s.low) .+ s.low
end

function contains(s::BoxSpace43{T,N}, v::AbstractVector{T}) where {T,N}
    all(v .>= s.low) && all(v .<= s.high)
end


################################################################################
# Pendulum
################################################################################

abstract type Environment43{T,M,A,O} end

mutable struct PendulumEnvironment43{T,M,A,O} <: Environment43{T,M,A,O}
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

function PendulumEnvironment43(; max_speed::T = 8.0, max_torque::T = 2.0,
        dt::T = 0.05, g::T = -10.0, m::T = 1.0, l::T = 1.0, s::Int = 1, vis::Visualizer = Visualizer()) where {T}
    mechanism = getmechanism(:pendulum, Δt = dt, g = g, m = m, l = l)
    nx = minCoordDim(mechanism)
    nu = controldim(mechanism)
    no = 3

    high = [1.0, 1.0, max_speed]
    aspace = BoxSpace27(controldim(mechanism), low = [-max_torque], high = [max_torque])
    ospace = BoxSpace27(no, low = -high, high = high)
    rng = MersenneTwister(s)
    x = Inf * ones(nx)
    last_u = Inf * ones(nu)
    build_robot!(mechanism, vis)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = PendulumEnvironment43{TYPES...}(mechanism, aspace, ospace, x, last_u,
        max_speed, nx, nu, no, max_torque, dt, g, m, l, rng, vis)
    seed(env, s = s)
    return env
end

function seed(env::PendulumEnvironment43{T}; s = 0) where {T}
    env.rng = MersenneTwister(s)
    return nothing
end

function reset(env::PendulumEnvironment43{T}) where {T}
    high = [π, 1.0]
    low = -high

    env.x = rand(env.rng, env.nx) .* (high .- low) .+ low
    env.last_u .= Inf
    return _get_obs(env)
end

function _get_obs(env::PendulumEnvironment43{T}) where {T}
    θ, ω = env.x
    return [cos(θ), sin(θ), ω]
end

function step(env::PendulumEnvironment43{T}, u::AbstractVector{T}) where {T}
    mechanism = env.mechanism
    x0 = env.x

    u0 = clamp.(u, -env.max_torque, env.max_torque)
    self.last_u = u0  # for rendering

    z0 = min2max(mechanism, x0)
    z1 = simon_step!(mechanism, z0, u0; ϵ = 1e-6, newtonIter = 100,
        lineIter = 10, verbose = false, btol = 1e-6, undercut = Inf)
    env.x = max2min(mechanism, z1)

    # Compute cost function
    θ0, ω0 = x0
    costs = angle_normalize(θ0)^2 + 1e-1 * ω0^2 + 1e-3 * u[1]^2 # angle_normalize enforces angle ∈ [-π, π]

    info = Dict()
    return _get_obs(env), -costs, False, info
end

function angle_normalize(x)
    return ((x + π) % (2 * π)) - π
end

function render(env::PendulumEnvironment43, mode="human"):
    z = min2max(env.mechanism, env.x)
    setrobot!(env.mechanism, env.vis, z)
    return nothing
end

function close(env::PendulumEnvironment43{M}; kwargs...) where {M}
    # visualizer stuff
    # if env.vis:
    #     env.vis.close()
    #     env.vis = None
    return nothing
end


#
# mech = getmechanism(:pendulum)
# minCoordDim(mech)
#
# env = make("Pendulum", dt = 0.01)
#
#
# reset(env)
# reset(env)
# reset(env)
# reset(env)
