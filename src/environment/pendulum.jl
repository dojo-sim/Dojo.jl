include("environment.jl")

################################################################################
# Pendulum
################################################################################

mutable struct PendulumEnvironment{T,M,A,O} <: Environment{T,M,A,O}
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

function PendulumEnvironment(; max_speed::T = 8.0, max_torque::T = 2.0,
        dt::T = 0.05, g::T = -10.0, m::T = 1.0, l::T = 1.0, s::Int = 1, vis::Visualizer = Visualizer()) where {T}
    mechanism = getmechanism(:pendulum, Δt = dt, g = g, m = m, l = l, damper = 0.5)
    nx = minCoordDim(mechanism)
    nu = controldim(mechanism)
    no = 3

    high = [1.0, 1.0, max_speed]
    aspace = BoxSpace(controldim(mechanism), low = [-max_torque], high = [max_torque])
    ospace = BoxSpace(no, low = -high, high = high)
    rng = MersenneTwister(s)
    x = Inf * ones(nx)
    last_u = Inf * ones(nu)
    # build_robot(vis, mechanism)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = PendulumEnvironment{TYPES...}(mechanism, aspace, ospace, x, last_u, nx, nu, no,
        max_speed, max_torque, dt, g, m, l, rng, vis)
    seed(env, s = s)
    return env
end

function seed(env::PendulumEnvironment{T}; s = 0) where {T}
    env.rng = MersenneTwister(s)
    return nothing
end

function reset(env::PendulumEnvironment{T}; x = nothing) where {T}
    if x != nothing
        env.x = x
    else
        high = [π, 1.0]
        low = -high
        @show env.nx
        @show size(rand(env.rng, env.nx))
        env.x = rand(env.rng, env.nx) .* (high .- low) .+ low
        env.last_u .= Inf
    end
    return _get_obs(env)
end

function _get_obs(env::PendulumEnvironment{T}) where {T}
    θ, ω = env.x
    return [cos(θ), sin(θ), ω]
end

function step(env::PendulumEnvironment{T}, u::AbstractVector{T}) where {T}
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
    costs = angle_normalize(θ0)^2 + 1e-1 * ω0^2 + 1e-3 * u[1]^2 # angle_normalize enforces angle ∈ [-π, π]

    info = Dict()
    return _get_obs(env), -costs, false, info
end

function angle_normalize(x)
    return ((x + π) % (2 * π)) - π
end

function render(env::PendulumEnvironment, mode="human")
    z = min2max(env.mechanism, env.x)
    set_robot(env.vis, env.mechanism, z)
    return nothing
end

function close(env::PendulumEnvironment{M}; kwargs...) where {M}
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



################################################################################
# Script
################################################################################

env = make("Pendulum", dt = 0.05, g = -9.81);
typeof(env)
reset(env, x = [0.1,0.2])
u = ones(1)
step(env, u)
minCoordDim(env.mechanism)

env = make("Pendulum");
for i = 1:20
    observation = reset(env)
    for t = 1:200
        # render(env)
        println(scn.(observation))
        u = sample(env.aspace)
        @show u
        @show env.x[2]
        observation, reward, done, info = step(env, u)
        if done
            println("Episode finished after $(t+1) timesteps")
            break
        end
    end
end
close(env)
