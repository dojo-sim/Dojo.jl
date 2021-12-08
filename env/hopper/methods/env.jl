################################################################################
# Hopper
################################################################################

struct Hopper{T,M,A,O} <: Environment{T,M,A,O} 
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

function Hopper(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer()) where T

    mechanism = gethopper(Δt=dt, g=g)
    initializehopper!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max 
        nx = maxCoordDim(mechanism)
    end
    nu = 3
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e-3 * ones(nu)), high=(1.0e-3 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))
    rng = MersenneTwister(s)
    z = hopper_nominal_max()
    x = mode == :min ? max2min(mechanism, z) : z
    fx = zeros(nx, nx) 
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = hopper_control_mask() 

    build_robot(vis, mechanism)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace)]
    env = Hopper{TYPES...}(mechanism, mode, aspace, ospace, 
        x, fx, fu,
        u_prev, control_mask,
        nx, nu, no,
        rng, vis)
    return env
end

function reset(env::Hopper; x=nothing)
    initializehopper!(env.mechanism) 

    if x != nothing
        env.x = x
    else
        z = hopper_nominal_max() 
        if env.mode == :min 
            env.x .= max2min(env.mechanism, z) 
        elseif env.mode == :max 
            env.x .= z 
        end
        env.u_prev .= 0.0
    end
    return _get_obs(env)
end

function _get_obs(env::Hopper)
    return env.x
end

function step(env::Hopper, x, u; diff=false)
    mechanism = env.mechanism
    Δt = mechanism.Δt

    x0 = x
    env.u_prev .= u  # for rendering

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, env.control_mask' * u; ϵ=1e-5, newtonIter=100,
        lineIter=10, verbose=false, btol=1e-5, undercut=1.5)
    env.x .= env.mode == :min ? max2min(mechanism, z1) : z1

    # Compute cost function
    costs = 0.0

    # Gradients
    if diff 
        if env.mode == :min 
            fx, fu = getMinGradients!(env.mechanism, z0, env.control_mask' * u, ϵ=1e-5, btol=1e-3, undercut=1.5, verbose=false)
        elseif env.mode == :max 
            fx, fu = getMaxGradients!(env.mechanism, z0, env.control_mask' * u, ϵ=1e-5, btol=1e-3, undercut=1.5, verbose=false)
        end
        env.fx .= fx
        env.fu .= fu * env.control_mask'
    end

    info = Dict()

    return _get_obs(env), -costs, false, info
end

function render(env::Hopper, mode="human")
    z = env.mode == :min ? min2max(env.mechanism, env.x) : env.x
    set_robot(env.vis, env.mechanism, z)
    return nothing
end

function close(env::Hopper; kwargs...) 
    return nothing
end

function hopper_nominal_max()
    # initial state
    x1b1 = [0.0; 0.0; 0.5]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; 0.0]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function hopper_control_mask()
    [0 0 0 1 0 0 0;
	 0 0 0 0 1 0 0;
	 0 0 0 0 0 0 1]
end