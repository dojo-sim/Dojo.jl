################################################################################
# Imports
################################################################################

import Base.contains
import Base.reset
import Base.step
import MeshCat.render

################################################################################
# Environment
################################################################################

type2symbol(H) = Symbol(lowercase(String(H.name.name)))

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
	control_mask::Matrix{T}
    control_scaling::Diagonal{T,Vector{T}} # could be merged with control mask
    nx::Int
    nu::Int
    no::Int
    info::I
    rng::Vector{MersenneTwister}
    vis::Visualizer
    opts_step::SolverOptions{T}
    opts_grad::SolverOptions{T}
end

function reset(env::Environment{X}; x=nothing) where X
    initialize!(env.mechanism, type2symbol(X))
    if x != nothing
        env.x = x
    else
        if env.mode == :min
            env.x .= get_minimal_state(env.mechanism)
        elseif env.mode == :max
            env.x .= get_maximal_state(env.mechanism)
        end
        env.u_prev .= 0.0
    end
    return _get_obs(env)
end

function _get_obs(env::Environment)
    return env.x
end

is_done(env::Environment, x) = false

function step(env::Environment, x, u; diff=false)
    mechanism = env.mechanism
    timestep = mechanism.timestep

    x0 = x
    # u = clip(env.aspace, u) # control limits
    env.u_prev .= u  # for rendering in Gym
	u_scaled = env.control_mask' * env.control_scaling * u

    z0 = env.mode == :min ? minimal_to_maximal(mechanism, x0) : x0
    z1 = step!(mechanism, z0, u_scaled; opts=env.opts_step)
    env.x .= env.mode == :min ? maximal_to_minimal(mechanism, z1) : z1

    # Compute cost
    costs = cost(env, x, u)

	# Check termination
	done = is_done(env, x)

    # Gradients
    if diff
        if env.mode == :min
            fx, fu = get_minimal_gradients(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        elseif env.mode == :max
            fx, fu = get_maximal_gradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        end
        env.fx .= fx
        env.fu .= fu * env.control_mask' * env.control_scaling
    end

    info = Dict()
    return _get_obs(env), -costs, done, info
end

function render(env::Environment, mode="human")
    z = env.mode == :min ? minimal_to_maximal(env.mechanism, env.x) : env.x
    set_robot(env.vis, env.mechanism, z, name=:robot)
    return nothing
end

function seed(env::Environment; s=0)
    env.rng[1] = MersenneTwister(s)
    return nothing
end

cost(env::Environment, x, u) = 0.0

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

function BoxSpace(n::Int; low::AbstractVector{T} = -ones(n), high::AbstractVector{T} = ones(n)) where T
    return BoxSpace{T,n}(n, low, high, (n,), T)
end

function sample(s::BoxSpace{T,N}) where {T,N}
    return rand(T,N) .* (s.high .- s.low) .+ s.low
end

function contains(s::BoxSpace{T,N}, v::AbstractVector{T}) where {T,N}
    all(v .>= s.low) && all(v .<= s.high)
end

function clip(s::BoxSpace, u)
    clamp.(u, s.low, s.high)
end

################################################################################
# Step
################################################################################
step(env::Environment, u; diff::Bool=false) = step(env, env.x, u; diff=diff)

function f(y, env::Environment, x, u, w)
	step(env, x, u)[1]
    y .= env.x
end

function fx(dx, env::Environment, x, u, w)
	step(env, x, u, diff=true)
    dx .= env.fx
end

function fu(du, env::Environment, x, u, w)
	# step(env, x, u, diff=true) # this is run in fx
	du .= env.fu
end

################################################################################
# Environments
# ##############################################################################
include("ant/methods/env.jl")
include("atlas/methods/env.jl")
include("cartpole/methods/env.jl")
include("halfcheetah/methods/env.jl")
include("hopper/methods/env.jl")
include("pendulum/methods/env.jl")
include("quadruped/methods/env.jl")
include("raiberthopper/methods/env.jl")
include("walker2d/methods/env.jl")
include("box/methods/env.jl")

################################################################################
# Visualize Trajectories
# ##############################################################################
function visualize(env::Environment, traj::Vector{Vector{T}}) where T
	@assert size(traj[1]) == size(env.x)
    storage = generate_storage(env.mechanism, [env.mode == :min ? minimal_to_maximal(env.mechanism, x) : x for x in traj])
    visualize(env.mechanism, storage, vis=env.vis)
end
