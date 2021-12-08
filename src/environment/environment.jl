
################################################################################
# Environment
################################################################################

mutable struct Environment32{M,A,O}
    mechanism::M
    model::Symbol
    reward_fct::Any # reward function R(s, a) = r
    iter::Int # current time step count
    max_iter::Int # number of time steps per episode
    aspace::A
    ospace::O
    # vis::Visualizer
end


function Environment32(mechanism::Mechanism, model::Symbol, reward_fct::Any; max_iter::Int = 100)
    iter = 0
    aspace = BoxSpace32(controldim(mechanism))
    ospace = BoxSpace32(maxCoordDim(mechanism))
    TYPE = [typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = Environment32{TYPE...}(mechanism, model, reward_fct, iter, max_iter, aspace, ospace)
    return env
end

function make(model::String; kwargs...)
    model = Symbol(model)
    mechanism = getmechanism(model; kwargs...)
    reward_fct = getreward(model)
    env = Environment32(mechanism, model, reward_fct)
    reset(env)
    return env
end

function reset(env::Environment32{M}; kwargs...) where {M}
    env.iter = 0
    initialize!(env.mechanism, Symbol(env.model); kwargs...)
    observation = getMaxState(env.mechanism)
    return observation
end

function render(env::Environment32{M}; kwargs...) where {M}

    return nothing
end

function close(env::Environment32{M}; kwargs...) where {M}
    # visualizer stuff
    return nothing
end

function need_to_reset(env::Environment32)
    done = env.iter >= env.max_iter - 1
    return done
end

function step(env::Environment32{M}, action) where {M}
    env.iter += 1
    mechanism = env.mechanism
    state = getMaxState(mechanism)
    observation = getMaxState(mechanism)

    reward = env.reward_fct(state, action)
    done = need_to_reset(env)
    info = Dict()
    return observation, reward, done, info
end



#
#
# ################################################################################
# # Script
# ################################################################################
#
# env = make("dice", Δt = 0.05, g = -9.81);
# reset(env, x = [0,0,0.2])
# action = nothing
# step(env, action)
#
#
# env = make("dice");
# for i = 1:20
#     observation = reset(env)
#     for t = 1:200
#         # render(env)
#         println(scn.(observation))
#         action = sample(env.aspace)
#         observation, reward, done, info = step(env, action)
#         if done
#             println("Episode finished after $(t+1) timesteps")
#             break
#         end
#     end
# end
# close(env)
#
# ################################################################################
# # Development
# ################################################################################
# # Utils
# function module_dir()
#     return joinpath(@__DIR__, "..", "..")
# end
#
# # Activate package
# using Pkg
# Pkg.activate(module_dir())
#
# # Load packages
# using Plots
# using Random
# using MeshCat
#
# # Open visualizer
# vis = Visualizer()
# open(vis)
#
# # Include new files
# include(joinpath(module_dir(), "examples", "loader.jl"))
#
# mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, conetype = :soc)
# # mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, conetype = :linear)
# # mech = getmechanism(:dice, Δt = 0.01, g = -9.81, contact = true, mode=:box, conetype = :impact)
# Random.seed!(100)
# ω = 0.0 * (rand(3) .- 0.5) * 1
# x = [0, 0, 1.0]
# v = 1.0 * [4, 1, 0.0]
# initialize!(mech, :dice, x = x, v = v, ω = ω)
# storage = simulate!(mech, 1.3, record = true, solver = :mehrotra!, verbose = false)
# visualize(mech, storage, vis = vis)
