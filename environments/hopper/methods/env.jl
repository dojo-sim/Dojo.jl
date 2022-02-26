################################################################################
# Hopper
################################################################################
struct Hopper end

function hopper(; 
    representation=:minimal, 
    timestep=0.05, 
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=1.0, 
    spring=10.0, 
    damper=50.0,
    seed=1, 
    contact_foot=true, 
    contact_body=true, 
    info=nothing, 
    vis=Visualizer(), 
    name=:robot,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5), 
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_hopper(
        timestep=timestep, 
        gravity=gravity, 
        friction_coefficient=friction_coefficient, 
        spring=spring, 
        damper=damper, 
        contact_foot=contact_foot,
        contact_body=contact_body)

    initialize_hopper!(mechanism)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu = 3
    no = nx

    # values taken from Mujoco's model, combining the control range -1, 1 and the motor gears.
    aspace = BoxSpace(nu, 
        low=(-Inf * ones(nu)), 
        high=(Inf * ones(nu)))
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
    motor_gear = [200, 200, 200.]
    control_scaling = Diagonal(timestep * motor_gear)

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [Hopper, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function reset(env::Environment{Hopper}; 
    x=nothing, reset_noise_scale=0.005)

    if x != nothing
        env.state .= x
    else
        # initialize above the ground to make sure that with random initialization we do not violate the ground constraint.
        initialize!(env.mechanism, :hopper, 
            z=0.25)
        x0 = get_minimal_state(env.mechanism)
        nx = minimal_dimension(env.mechanism)

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

function get_observation(env::Environment{Hopper}; 
    full_state=false)

    full_state && (return env.state)

    nx = minimal_dimension(env.mechanism)
    if env.representation == :minimal
        o = env.state
    elseif env.representation == :maximal
        o = maximal_to_minimal(env.mechanism, env.state)
    end

    # clamp velocities and remove x position
    ind = velocity_index(env.mechanism)
    vel = clamp.(o[ind], -10, 10)
    o[ind] .= vel
    o = o[[1; 3:nx]]
    return o
end

function cost(env::Environment{Hopper}, x, u;
        alive_bonus=0.1)

    if env.representation == :minimal
        x_velocity = -x[5]
    else
        i_torso = findfirst(body -> body.name == "torso", collect(env.mechanism.bodies))
        z_torso = x[(i_torso-1)*13 .+ (1:13)]
        x_velocity = -z_torso[4]
    end
    # c = -x_velocity/10 # Gym reward
    c = -x_velocity
    c -= alive_bonus
    c += 1e-4 * u'*u
    return c
end

function is_done(::Environment{Hopper}, x)
    nx = minimal_dimension(env.mechanism)
    if env.representation == :minimal
        x0 = x
    elseif env.representation == :maximal
        x0 = maximal_to_minimal(env.mechanism, x)
    end
    height = x0[1]
    ang = x0[3]
    done = !(
        all(isfinite.(x0)) &&
        all(abs.(x0[[1;3:nx]]) .< 100) &&
        (height > 0.7) &&
        (abs(ang) < 0.2)
        )
    return done
end

