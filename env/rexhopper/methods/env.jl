################################################################################
# Hopper
################################################################################
struct RexHopper end

function rexhopper(;
    mode::Symbol=:min,
    model=:rexhopper,
    timestep::T=0.05,
    s=1,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient::T=1.0,
    spring=0.0,
    damper=1.0,
    contact_type=:nonlinear,
    contact::Bool=true,
    contact_body::Bool=true,
    limits::Bool = true,
    info=nothing,
    infeasible_control::Bool=false,
    vis::Visualizer=Visualizer(),
    name::Symbol=:robot,
    opts_step=SolverOptions(rtol=1.0e-4, btol=1.0e-4, undercut=10.0),
    opts_grad=SolverOptions(rtol=1.0e-4, btol=1.0e-4, undercut=10.0)) where T

    mechanism = get_rexhopper(model=model,
        limits=limits,
        contact_type=contact_type,
        timestep=timestep,
        gravity=gravity,
        friction_coefficient=friction_coefficient,
        spring=spring,
        damper=damper,
        contact=contact,
        contact_body=contact_body)

    initialize_rexhopper!(mechanism)

    if mode == :min
        nx = minimal_dimension(mechanism)
    elseif mode == :max
        nx = maximal_dimension(mechanism)
    end
    nu_inf = input_dimension(mechanism)
    nu = infeasible_control ? nu_inf : 5
    no = nx

    # values taken from Mujoco's model, combining the control range -1, 1 and the motor gears.
    aspace = BoxSpace(nu, low=(-Inf * ones(nu)), high=(Inf * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = get_maximal_state(mechanism)
    x = mode == :min ? maximal_to_minimal(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = infeasible_control ? I(nu) : cat(zeros(3,3), I(3), 1, 0, 1, zeros(5,5), dims=(1,2))
    motor_gear = ones(nu)
    control_scaling = Diagonal(motor_gear)

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [RexHopper, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function reset(env::Environment{RexHopper}; x=nothing)
    if x != nothing
        env.x .= x
    else
        # initialize above the ground to make sure that with random initialization we do not violate the ground constraint.
        initialize!(env.mechanism, :rexhopper)
        x0 = get_minimal_state(env.mechanism)
        nx = minimal_dimension(env.mechanism)
        z = minimal_to_maximal(env.mechanism, x)
        set_state!(env.mechanism, z)
        if env.mode == :min
            env.x .= get_minimal_state(env.mechanism)
        elseif env.mode == :max
            env.x .= get_maximal_state(env.mechanism)
        end
        env.u_prev .= 0.0
    end
    return _get_obs(env)
end
