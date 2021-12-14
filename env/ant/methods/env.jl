################################################################################
# Ant
################################################################################
struct Ant end

function ant(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    cf::T=0.5, spring::T=50.0, damper::T=50.0, s::Int=1,
    contact::Bool=true, contact_body=true,
    info=nothing, vis::Visualizer=Visualizer(),
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = getant(Δt=dt, g=g, cf=cf, spring=spring, damper=damper,contact=contact, contact_body=contact_body)
    initializeant!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max
        nx = maxCoordDim(mechanism)
    end
    nu = 8
    no = nx

    aspace = BoxSpace(nu, low=(-ones(nu)), high=(ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(8, 6) I(nu)]
    control_scaling = Diagonal(dt * 150.0 * ones(nu))

    build_robot(vis, mechanism)

    TYPES = [Ant, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
    return env
end

function step(env::Environment{Ant}, x, u; diff=false)
    mechanism = env.mechanism
    Δt = mechanism.Δt

    # x position (before)
    xposbefore = x[1]

    x0 = copy(x)
    env.u_prev .= u  # for rendering in Gym
	u_scaled = env.control_mask' * env.control_scaling * u

    z0 = env.mode == :min ? min2max(mechanism, x0) : x0
    z1 = step!(mechanism, z0, u_scaled; opts=env.opts_step)
    env.x .= env.mode == :min ? max2min(mechanism, z1) : z1

    # x position (after)
    xposafter = env.x[1]

    # forward reward
    forward_reward = 2*(xposafter - xposbefore) / Δt

    # control cost
    # ctrl_cost = (0.5 * u' * u)[1]
	ctrl_cost = (0.05 * u' * u)[1]

    # contact cost
    contact_cost = 0.0

    for ineq in mechanism.ineqconstraints
        contact_cost += 0.5 * 1.0e-3 * max(-1.0, min(1.0, ineq.γsol[2][1]))^2.0
    end

    # survive reward
	# survive_reward = 1.0
    survive_reward = 0.05

    # total reward
    reward = forward_reward - ctrl_cost - contact_cost + survive_reward

    # done ?
    done = !(all(isfinite.(env.x)) && (env.x[3] >= 0.2) && (env.x[3] <= 1.0))

    # Gradients
    if diff
        if env.mode == :min
            fx, fu = getMinGradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        elseif env.mode == :max
            fx, fu = getMaxGradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        end
        env.fx .= fx
        env.fu .= fu * env.control_mask' * env.control_scaling
    end

    info = Dict()

    return _get_obs(env), reward, done, info
end

function reset(env::Environment{Ant};
    x=nothing,
    pos_noise=Distributions.Uniform(-0.1, 0.1),
    vel_noise=Distributions.Normal(0.0, 0.01))

    initialize!(env.mechanism, type2symbol(Ant))

    if x != nothing
        env.x .= x
    else
        x = getMinState(env.mechanism, pos_noise=pos_noise, vel_noise=vel_noise)
        if env.mode == :min
            setState!(env.mechanism, min2max(env.mechanism, x))
            env.x .= x
        elseif env.mode == :max
            z = min2max(env.mechanism, x)
            setState!(env.mechanism, z)
            env.x .= z
        end
        env.u_prev .= 0.0
    end

    return _get_obs(env)
end

function _get_obs(env::Environment{Ant,T}) where T
    contact_force = T[]
    for ineq in env.mechanism.ineqconstraints
        push!(contact_force, max(-1.0, min(1.0, ineq.γsol[2][1])))
    end
    return [env.x; contact_force]
end

# env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=25.0, spring=25.0, contact=true, contact_body=true)
# total_mass(env.mechanism)

# # initialize!(env.mechanism, :ant)
# open(env.vis)
# # storage = simulate!(env.mechanism, 1.0, record=true, verbose=false)
# # visualize(env.mechanism, storage, vis=env.vis)

# reset(env)
# render(env)
# x0 = getMinState(env.mechanism)

# for i = 1:100
#     u = 1.0 * rand(Distributions.Uniform(-1.0, 1.0), 8)
#     x0, r, _ = step(env, x0, u)
#     @show r
#     render(env)
# end
