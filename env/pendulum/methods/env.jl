################################################################################
# Pendulum
################################################################################
struct Pendulum end

function pendulum(; mode::Symbol=:min, max_speed::T=8.0, max_torque::T=8.0,
        dt::T=0.05, g::T=-10.0, m::T=1.0, l::T=1.0, damper=0.0, s::Int=1, vis::Visualizer=Visualizer(),
        opts_step::InteriorPointOptions = InteriorPointOptions(),
        opts_grad::InteriorPointOptions = InteriorPointOptions()) where {T}

    mechanism = getmechanism(:pendulum, Δt=dt, g=g, m=m, l=l, damper=damper)
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
    rng = [MersenneTwister(s),]

    x = Inf * ones(nx)
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = Inf * ones(nu)
    control_mask = ones(1,1)
    build_robot(vis, mechanism)

    info = Dict(:max_speed => max_speed, :max_torque => max_torque)

    TYPES = [T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{Pendulum, TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, nx, nu, no, info,
        rng, vis, opts_step, opts_grad)
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
            env.x .= rand(env.rng[1], env.nx) .* (high .- low) .+ low
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
    max_torque = env.info[:max_torque]

    x0 = x
    u0 = clamp.(u, -max_torque, max_torque)
    env.u_prev .= u0  # for rendering

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, Δt * u0; opts = env.opts_step)
    env.x .= env.mode == :min ? max2min(mechanism, z1) : z1

    # Compute cost function
    costs = cost(env, x0, u0)

    # Gradients
    if diff
        if env.mode == :min
            fx, fu = getMinGradients!(env.mechanism, z0, Δt * u0, opts=env.opts_grad)
        elseif env.mode == :max
            fx, fu = getMaxGradients!(env.mechanism, z0, Δt * u0, opts=env.opts_grad)
        end
        env.fx .= fx
        env.fu .= Δt * fu
    end

    info = Dict()
    return _get_obs(env), -costs, false, info
end

function angle_normalize(x)
    return ((x + π) % (2 * π)) - π
end

function pendulum_nominal_max()
    x1 = [0.0; 0.0; -0.5]
    v1 = [0.0; 0.0; 0.0]
    q1 = [1.0; 0.0; 0.0; 0.0]
    ω1 = [0.0; 0.0; 0.0]
    z1 = [x1; v1; q1; ω1]
end

function pendulum_goal_max()
    xT = [0.0; 0.0; 0.5]
    vT = [0.0; 0.0; 0.0]
    qT = [0.0; 1.0; 0.0; 0.0]
    ωT = [0.0; 0.0; 0.0]
    zT = [xT; vT; qT; ωT]
end

function cost(env::Environment{Pendulum}, x, u)
    if env.mode == :min
        θ, ω = x
        c = angle_normalize(θ)^2 + 1e-1 * ω^2 + 1e-3 * u[1]^2 # angle_normalize enforces angle ∈ [-π, π]
    else
        c = Inf
    end
    return c
end
