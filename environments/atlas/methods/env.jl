"""
    Atlas <: Environment

    humanoid robot designed by Boston Dynamics (v5)
"""
struct Atlas end

function atlas(; 
    representation=:minimal, 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=0.8,
    damper=10.0, 
    spring=0.0, 
    info=nothing, 
    model_type=:simple,
    seed=1, 
    contact_feet=true, 
    contact_body=false,
    vis=Visualizer(), 
    name=:robot,
    infeasible_control=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5)) where T

    mechanism = get_mechanism(:atlas, 
        timestep=timestep, 
        gravity=gravity, 
        friction_coefficient=friction_coefficient, 
        damper=damper,
        spring=spring, 
        contact_feet=contact_feet,
        contact_body=contact_body,
        model_type=model_type)

    initialize!(mechanism, :atlas)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu_inf = input_dimension(mechanism)
    nu = infeasible_control ? nu_inf : nu_inf - 6 # remove first 6 controls
    no = nx

    aspace = BoxSpace(nu, 
        low=(-1.0e8 * ones(nu)), 
        high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, 
        low=(-Inf * ones(no)), 
        high=(Inf * ones(no)))

    rng = MersenneTwister(seed)
    z = get_maximal_state(mechanism)
    x = representation == :minimal ? maximal_to_minimal(mechanism, z) : z
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = infeasible_control ? I(nu) : [zeros(nu, 6) I(nu)]
    control_scaling = Diagonal(ones(nu))

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [Atlas, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end
