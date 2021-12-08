################################################################################
# Hopper
################################################################################

struct Hopper{T,M,A,O} <: Environment{T,M,A,O} 
    mechanism::M
    mode::Symbol
    aspace::A
    ospace::O
    x::Vector{T}
    u_prev::Vector{T}
    nx::Int
    nu::Int
    no::Int
    rng::MersenneTwister
    vis::Visualizer
end

function Hopper(; mode::Symbol=:min, dt::T=0.05, g::T=-10.0, m::T=1.0, l::T=1.0, 
    damper=0.5, s::Int=1, vis::Visualizer=Visualizer()) where T
    mechanism = getmechanism(:pendulum, Δt=dt, g=g, m=m, l=l, damper=damper)
    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max 
        nx = maxCoordDim(mechanism)
    end
    nu = controldim(mechanism)
    no = nx

    aspace = BoxSpace(controldim(mechanism), low=[-Inf, -Inf, -Inf], high=[Inf, Inf, Inf])
    ospace = BoxSpace(no, low=-high, high=high)
    rng = MersenneTwister(s)
    x = Inf * ones(nx)
    u_prev = Inf * ones(nu)
    build_robot(vis, mechanism)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = Hopper{TYPES...}(mechanism, mode, aspace, ospace, x, u_prev,
        nx, nu, no,
        rng, vis)
    seed(env, s=s)
    return env
end

function seed(env::Hopper; s=0)
    env.rng = MersenneTwister(s)
    return nothing
end

function reset(env::Hopper; x=nothing)
    if x != nothing
        env.x = x
    else
        high = [π, 1.0]
        low = -high
        env.x = rand(env.rng, env.nx) .* (high .- low) .+ low
        env.u_prev .= Inf
    end
    return _get_obs(env)
end

function _get_obs(env::Hopper)
    return env.x
end

function step(env::Hopper, u)
    mechanism = env.mechanism
    Δt = mechanism.Δt
    x0 = env.x

    env.u_prev = u  # for rendering

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, u0; ϵ=1e-6, newtonIter=100,
        lineIter=10, verbose=false, btol=1e-6, undercut=Inf)
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

function render(env::Hopper, mode="human")
    z = min2max(env.mechanism, env.x)
    set_robot(env.vis, env.mechanism, z)
    return nothing
end

function close(env::Hopper; kwargs...) 
    # visualizer stuff
    # if env.vis:
    #     env.vis.close()
    #     env.vis = None
    return nothing
end