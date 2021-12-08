################################################################################
# Pendulum
################################################################################
struct Pendulum end 

function pendulum(; mode::Symbol=:min, max_speed::T=8.0, max_torque::T=8.0,
        dt::T=0.05, g::T=-10.0, m::T=1.0, l::T=1.0, s::Int=1, vis::Visualizer=Visualizer(),
        info=nothing,
        opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where {T}

    mechanism = getmechanism(:pendulum, Δt=dt, g=g, m=m, l=l, damper=0.5)
    initialize!(mechanism, :pendulum)

    if mode == :min
        nx = minCoordDim(mechanism)
        no = 3
    elseif mode == :max 
        nx = maxCoordDim(mechanism) 
        no = 13
    end
    nu = controldim(mechanism)

    high = [1.0, 1.0, max_speed]
    aspace = BoxSpace(controldim(mechanism), low=[-max_torque], high=[max_torque])
    ospace = BoxSpace(no, low=-high, high=high)
    rng = MersenneTwister(s)

    x = Inf * ones(nx)
    fx = zeros(nx, nx) 
    fu = zeros(nx, nu) 

    u_prev = Inf * ones(nu)
    build_robot(vis, mechanism)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{Pendulum, TYPES...}(mechanism, mode, aspace, ospace, 
        x, fx, fu,
        u_prev, nx, nu, no,
        rng, vis,
        info,
        opts_step, opts_grad)

    return env
end

function reset(env::Environment{Pendulum}; x=nothing)
    initialize!(env.mechanism, :pendulum)

    if x != nothing
        env.x .= x
    else
        if env.mode == :min 
            high = [π, 1.0]
            low = -high
            env.x .= rand(env.rng, env.nx) .* (high .- low) .+ low
        elseif env.mode == :max 
            env.x .= pendulum_nominal_max()
        end
        env.u_prev .= Inf
    end
    return _get_obs(env)
end

function _get_obs(env::Environment{Pendulum})
    if env.mode == :min
        θ, ω = env.x
        return [cos(θ), sin(θ), ω]
    else env.mode == :max 
        return env.x 
    end
end

function step(env::Environment{Pendulum}, x, u; diff=false)
    mechanism = env.mechanism
    Δt = mechanism.Δt

    x0 = x
    u0 = clamp.(u, -env.max_torque, env.max_torque)
    env.u_prev .= u0  # for rendering

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, Δt * u0; opts=env.opts_step)
    env.x .= env.mode == :min ? max2min(mechanism, z1) : z1

    # Compute costs
    costs = reward(env, x0, u0)

    # Gradients
    if diff 
        if env.mode == :min 
            fx, fu = getMinGradients!(env.mechanism, z0, u, opts=env.opts_grad)
        elseif env.mode == :max 
            fx, fu = getMaxGradients!(env.mechanism, z0, u, opts=env.opts_grad)
        end
        env.fx .= fx
        env.fu .= fu 
    end

    info = Dict()
    return _get_obs(env), -costs, false, info
end

function angle_normalize(x)
    return ((x + π) % (2 * π)) - π
end

function reward(env::Environment{Pendulum}, x, u)
    x0 = env.mode == :max ? max2min(env.mechanism, x) : x 
    θ0, ω0 = x0
    costs = angle_normalize(θ0)^2 + 1e-1 * ω0^2 + 1e-3 * u[1]^2 # angle_normalize enforces angle ∈ [-π, π]
    return costs 
end