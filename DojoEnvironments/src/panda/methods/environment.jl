"""
    Panda <: Environment

    6-degree-of-freedom robotic arm designed by Franka Emika
"""
struct Panda end

function panda(;
    representation=:minimal,
    timestep=0.01,
    gravity=-9.81,
    friction_coefficient=0.8,
    dampers=0,
    springs=0,
    parse_dampers=true,
    info=nothing,
    model_type=:end_effector,
    seed=1,
    contact=true,
    limits=true,
    vis=Visualizer(),
    name=:robot,
    # infeasible_control=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_mechanism(:panda;
        timestep,
        gravity,
        friction_coefficient,
        damper,
        spring,
        parse_damper,
        contact,
        limits,
        model_type)

    initialize!(mechanism, :panda)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu = input_dimension(mechanism)
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
    control_mask = I(nu)
    control_scaling = Diagonal(ones(nu))

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [Panda, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end

# function initialize_panda!(mechanism::Mechanism{T};
#     joint_angles=[[0,-0.8,0,1.6,0,-2.4,0]; zeros(input_dimension(mechanism)-7)],
#     joint_velocities=zeros(input_dimension(mechanism))) where T

#     nu = input_dimension(mechanism)
#     zero_velocity!(mechanism)
#     y = vcat([[joint_angles[i], joint_velocities[i]] for i=1:nu]...)
#     set_minimal_state!(mechanism, y)
#     return nothing
# end