################################################################################
# HalfCheetah
################################################################################
struct HalfCheetah end

function halfcheetah(; 
    representation=:minimal, 
    timestep=0.05, 
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.4, 
    spring=[240.0, 180.0, 120.0, 180.0, 120.0, 60.0], 
    damper=2.5 * [6., 4.5, 3., 4.5, 3., 1.5],
    limits=true,
    seed=1, 
    contact_feet=true, 
    contact_body=true,
    info=nothing, 
    vis=Visualizer(), 
    name=:robot,
    opts_step=SolverOptions(), 
    opts_grad=SolverOptions(),
    T=Float64)

    mechanism = get_halfcheetah(
        timestep=timestep, 
        gravity=gravity, 
        friction_coefficient=friction_coefficient, 
        spring=spring, 
        damper=damper, 
        contact_feet=contact_feet, 
        contact_body=contact_body, 
        limits=limits)
    initialize_halfcheetah!(mechanism)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu = 6
    no = nx

    # values taken from Mujoco's model, combining the control range -1, 1 and the motor gears.
    aspace = BoxSpace(nu, 
        low=(-ones(nu)), 
        high=(ones(nu)))
    ospace = BoxSpace(no, 
        low=(-Inf * ones(no)), 
        high=(Inf * ones(no)))

    rng = MersenneTwister(seed)

    z = get_maximal_state(mechanism)
    x = representation == :minimal ? maximal_to_minimal(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(nu, 3) I(nu)]
    motor_gear = [120, 90, 60, 120, 60, 30.]
    control_scaling = Diagonal(timestep * motor_gear)

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [HalfCheetah, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function reset(env::Environment{HalfCheetah}; 
    x=nothing, reset_noise_scale=0.1)
    if x != nothing
        env.state .= x
    else
        # initialize above the ground to make sure that with random initialization we do not violate the ground constraint.
        initialize!(env.mechanism, :halfcheetah, 
            z=0.25)
        x0 = get_minimal_state(env.mechanism)
        nx = minimal_dimension(env.mechanism)
        nz = maximal_dimension(env.mechanism)

        low = -reset_noise_scale
        high = reset_noise_scale
        x = x0 + (high - low) .* rand(env.rng[1], nx) .+ low # we ignored the normal distribution on the velocities
        z = minimal_to_maximal(env.mechanism, x)
        set_state!(env.mechanism, z)
        if env.representation == :minimal
            env.state .= get_minimal_state(env.mechanism)
        elseif env.representation == :maximal
            env.state .= get_maximal_state(env.mechanism)
        end
        env.input_previous .= 0.0
    end
    return get_observation(env)
end

function cost(env::Environment{HalfCheetah}, x, u;
        forward_reward_weight=1.0, 
        ctrl_cost_weight=0.1)

    if env.representation == :minimal
        x_velocity = -x[5]
    else
        i_torso = findfirst(body -> body.name == "torso", env.mechanism.bodies)
        z_torso = x[(i_torso-1) * 13 .+ (1:13)]
        x_velocity = -z_torso[4]
    end
    # @show mean(abs.(u))
    c = ctrl_cost_weight * u' * u - x_velocity * forward_reward_weight
    return c
end


