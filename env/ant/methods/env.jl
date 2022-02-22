################################################################################
# Ant
################################################################################
struct Ant end

function ant(; mode::Symbol=:min, dt::T=0.05, gravity=[0.0; 0.0; -9.81],
    friction_coefficient::T=0.5, spring=0.0, damper::T=1.0, s::Int=1,
    contact::Bool=true, contact_body=true,
    limits::Bool=true,
    info=nothing, vis::Visualizer=Visualizer(), name::Symbol=:robot,
    opts_step=SolverOptions(), opts_grad=SolverOptions()) where T

    mechanism = get_ant(timestep=dt, gravity=gravity, friction_coefficient=friction_coefficient, spring=spring, damper=damper,contact=contact, contact_body=contact_body, limits=limits)
    initialize_ant!(mechanism)

    if mode == :min
        nx = minimal_dimension(mechanism)
    elseif mode == :max
        nx = maximal_dimension(mechanism)
    end
    nu = 8
    no = nx

    aspace = BoxSpace(nu, low=(-ones(nu)), high=(ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = get_maximal_state(mechanism)
    x = mode == :min ? maximal_to_minimal(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(8, 6) I(nu)]
    control_scaling = Diagonal(dt * 150.0 * ones(nu))

    build_robot(mechanism, vis=vis, name=name)

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
    timestep = mechanism.timestep

    # x position (before)
    xposbefore = x[1]

    x0 = copy(x)
    env.u_prev .= u  # for rendering in Gym
	u_scaled = env.control_mask' * env.control_scaling * u

    z0 = env.mode == :min ? minimal_to_maximal(mechanism, x0) : x0
    z1 = step!(mechanism, z0, u_scaled; opts=env.opts_step)
    env.x .= env.mode == :min ? maximal_to_minimal(mechanism, z1) : z1

    # x position (after)
    xposafter = env.x[1]

    # forward reward
    forward_reward = 2*(xposafter - xposbefore) / timestep

    # control cost
    # ctrl_cost = (0.5 * u' * u)[1]
	ctrl_cost = (0.05 * u' * u)[1]

    # contact cost
    contact_cost = 0.0

    for contact in mechanism.contacts
        contact_cost += 0.5 * 1.0e-3 * max(-1.0, min(1.0, contact.impulses[2][1]))^2.0
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
            fx, fu = get_minimal_gradients(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        elseif env.mode == :max
            fx, fu = get_maximal_gradients!(env.mechanism, z0, u_scaled, opts=env.opts_grad)
        end
        env.fx .= fx
        env.fu .= fu * env.control_mask' * env.control_scaling
    end

    info = Dict()

    return _get_obs(env), reward, done, info
end

function reset(env::Environment{Ant};
    x=nothing,
    pos_noise=Uniform(-0.1, 0.1),
    vel_noise=Normal(0.0, 0.01))

    initialize!(env.mechanism, type2symbol(Ant))

    if x != nothing
        env.x .= x
    else
        x = get_minimal_state(env.mechanism, pos_noise=pos_noise, vel_noise=vel_noise)
        if env.mode == :min
            set_state!(env.mechanism, minimal_to_maximal(env.mechanism, x))
            env.x .= x
        elseif env.mode == :max
            z = minimal_to_maximal(env.mechanism, x)
            set_state!(env.mechanism, z)
            env.x .= z
        end
        env.u_prev .= 0.0
    end

    return _get_obs(env)
end

function _get_obs(env::Environment{Ant,T}) where T
    contact_force = T[]
    for contact in env.mechanism.contacts
        push!(contact_force, max(-1.0, min(1.0, contact.impulses[2][1])))
    end
    return [env.x; contact_force]
end

# env = make("ant", mode=:min, g=-9.81, dt=0.05, damper=1.0, spring=0.0, limits=true, contact=true, contact_body=true)
# initialize!(env.mechanism, :ant)
# open(env.vis)
# storage = simulate!(env.mechanism, 1.0, record=true, verbose=false)
# visualize(env.mechanism, storage, vis=env.vis, show_contact=true)

# reset(env)
# render(env)
# x0 = get_minimal_state(env.mechanism)

# for i = 1:100
#     u = 0.0 * rand(Distributions.Uniform(-1.0, 1.0), 8)
#     x0, r, _ = step(env, x0, u)
#     @show r
#     render(env)
# end
