################################################################################
# Block
################################################################################
struct Block end

function block(; mode::Symbol=:max, dt::T=0.05, gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8, side=0.5, contact=true, contact_type=:contact,
    s::Int=1, vis::Visualizer=Visualizer(), info=nothing, name::Symbol=:robot,
    control_scaling=Diagonal(ones(3)),
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5), opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5)) where T

    mechanism = get_box(timestep=dt, gravity=gravity, friction_coefficient=friction_coefficient, side=side, contact=contact, contact_type=contact_type)
    initialize_box!(mechanism)

    if mode == :min
        nx = minimal_dimension(mechanism)
    elseif mode == :max
        nx = maximal_dimension(mechanism)
    end
    nu = 3
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e8 * ones(nu)), high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = get_maximal_state(mechanism)
    x = mode == :min ? maximal_to_minimal(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [I(nu) zeros(nu, nu)]

    build_robot(vis, mechanism, name=name)

    TYPES = [Block, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end