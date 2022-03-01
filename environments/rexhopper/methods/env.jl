"""
    RexHopper <: Environment

    hopping robot designed and build by the Robotic Exploration Lab (Benjamin Bokser)
"""
struct RexHopper end

function rexhopper(;
    representation=:minimal,
    model=:rexhopper,
    timestep=0.05,
    seed=1,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=1.0,
    spring=0.0,
    damper=1.0,
    contact_type=:nonlinear,
    contact_foot=true,
    contact_body=true,
    limits=true,
    info=nothing,
    infeasible_control=false,
    vis=Visualizer(),
    name=:robot,
    opts_step=SolverOptions(rtol=1.0e-4, btol=1.0e-4, undercut=10.0),
    opts_grad=SolverOptions(rtol=1.0e-4, btol=1.0e-4, undercut=10.0),
    T=Float64)

    mechanism = get_rexhopper(
        model=model,
        limits=limits,
        contact_type=contact_type,
        timestep=timestep,
        gravity=gravity,
        friction_coefficient=friction_coefficient,
        spring=spring,
        damper=damper,
        contact_foot=contact_foot,
        contact_body=contact_body)

    initialize_rexhopper!(mechanism)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu_inf = input_dimension(mechanism)
    nu = infeasible_control ? nu_inf : 5
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
    control_map = infeasible_control ? 1.0 * I(nu) : cat(zeros(3, 3), 1.0 * I(3), 1.0, 0.0, 1.0, zeros(5, 5), dims=(1,2))

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [RexHopper, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, 
        representation, 
        aspace, 
        ospace,
        x, fx, fu,
        u_prev, control_map,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function reset(env::Environment{RexHopper}; 
    x=nothing)

    if x != nothing
        env.state .= x
    else
        # initialize above the ground to make sure that with random initialization we do not violate the ground constraint.
        initialize!(env.mechanism, :rexhopper)
        x0 = get_minimal_state(env.mechanism)
        nx = minimal_dimension(env.mechanism)
        z = minimal_to_maximal(env.mechanism, x)
        set_maximal_state!(env.mechanism, z)
        if env.representation == :minimal
            env.state .= get_minimal_state(env.mechanism)
        elseif env.representation == :maximal
            env.state .= get_maximal_state(env.mechanism)
        end
        env.input_previous .= 0.0
    end
    return get_observation(env)
end
