################################################################################
# Quadruped
################################################################################

struct Quadruped{T,M,A,O} <: Environment{T,M,A,O} 
    mechanism::M
    mode::Symbol
    aspace::A
    ospace::O
    x::Vector{T}
    fx::Matrix{T} 
    fu::Matrix{T}
    u_prev::Vector{T}
    control_mask::Matrix{T}
    nx::Int
    nu::Int
    no::Int
    rng::MersenneTwister
    vis::Visualizer
end

function Quadruped(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81, cf=0.8,
    damper=10.0, spring=0.0,
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer()) where T
    mechanism = getquadruped(Δt=dt, g=g, cf=cf, damper=damper, spring=spring)
    initialize!(mechanism, :quadruped)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max 
        nx = maxCoordDim(mechanism)
    end
    nu = 12
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e-3 * ones(nu)), high=(1.0e-3 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))
    rng = MersenneTwister(s)
    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z
    fx = zeros(nx, nx) 
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = quadruped_control_mask() 

    build_robot(vis, mechanism)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = Quadruped{TYPES...}(mechanism, mode, aspace, ospace, 
        x, fx, fu,
        u_prev, control_mask,
        nx, nu, no,
        rng, vis)
    return env
end

function reset(env::Quadruped; x=nothing)
    initialize!(env.mechanism, :quadruped) 

    if x != nothing
        env.x = x
    else
        z = getStateMax(env.mechanism) 
        if env.mode == :min 
            env.x .= max2min(env.mechanism, z) 
        elseif env.mode == :max 
            env.x .= z 
        end
        env.u_prev .= 0.0
    end
    return _get_obs(env)
end

function _get_obs(env::Quadruped)
    return env.x
end

function step(env::Quadruped, x, u; diff=false)
    mechanism = env.mechanism
    Δt = mechanism.Δt

    x0 = x
    env.u_prev .= u  # for rendering

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, env.control_mask' * u; ϵ=3e-4, newtonIter=100,
        lineIter=10, verbose=false, btol=3e-4, undercut=1.5)
    env.x .= env.mode == :min ? max2min(mechanism, z1) : z1

    # Compute cost function
    costs = 0.0

    # Gradients
    if diff 
        if env.mode == :min 
            fx, fu = getMinGradients!(env.mechanism, z0, env.control_mask' * u, ϵ=3e-4, btol=3e-4, undercut=1.5, verbose=false)
        elseif env.mode == :max 
            fx, fu = getMaxGradients!(env.mechanism, z0, env.control_mask' * u, ϵ=3e-4, btol=3e-4, undercut=1.5, verbose=false)
        end
        env.fx .= fx
        env.fu .= fu * env.control_mask'
    end

    info = Dict()

    return _get_obs(env), -costs, false, info
end

function render(env::Quadruped, mode="human")
    z = env.mode == :min ? min2max(env.mechanism, env.x) : env.x
    set_robot(env.vis, env.mechanism, z)
    return nothing
end

function close(env::Quadruped; kwargs...) 
    return nothing
end

function quadruped_control_mask()
    [zeros(12, 6) I(12)]
end