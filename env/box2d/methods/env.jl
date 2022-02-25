################################################################################
# Quadruped
################################################################################
struct Box2d end

function box2d(; mode::Symbol=:min, dt::T=0.05, gravity=[0.0; 0.0; -9.81], friction_coefficient=0.8,
    info=nothing,
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer(), name::Symbol=:robot,
    infeasible_control::Bool=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5)) where T

    mechanism = get_mechanism(:box2d, timestep=dt, gravity=gravity,
        friction_coefficient=friction_coefficient)
    initialize!(mechanism, :box2d)

    if mode == :min
        nx = minimal_dimension(mechanism)
    elseif mode == :max
        nx = maximal_dimension(mechanism)
    end
    nu_inf = input_dimension(mechanism)
    nu = infeasible_control ? nu_inf : nu_inf - 2 # remove last 2 controls
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e8 * ones(nu)), high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)
    z = get_maximal_state(mechanism)
    x = mode == :min ? maximal_to_minimal(mechanism, z) : z
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = infeasible_control ? I(nu) : [I(nu) zeros(nu, 2)]
    control_scaling = Diagonal(ones(nu))

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [Quadruped, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end
