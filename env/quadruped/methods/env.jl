################################################################################
# Quadruped
################################################################################
struct Quadruped end

function quadruped(; mode::Symbol=:min, dt::T=0.05, gravity=[0.0; 0.0; -9.81], friction_coefficient=0.8,
    damper=10.0, spring=0.0, info=nothing,
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer(), name::Symbol=:robot,
    infeasible_control::Bool=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5)) where T

    mechanism = get_mechanism(:quadruped, timestep=dt, gravity=gravity, friction_coefficient=friction_coefficient, damper=damper, spring=spring)
    initialize!(mechanism, :quadruped)

    if mode == :min
        nx = minimal_dimension(mechanism)
    elseif mode == :max
        nx = maximal_dimension(mechanism)
    end
    nu_inf = control_dimension(mechanism)
    nu = infeasible_control ? nu_inf : nu_inf - 6 # remove first 6 controls
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e8 * ones(nu)), high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)
    z = get_maximal_state(mechanism)
    x = mode == :min ? maximal_to_minimal(mechanism, z) : z
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = infeasible_control ? I(nu) : [zeros(nu, 6) I(nu)]
    control_scaling = Diagonal(ones(nu))

    build_robot(vis, mechanism, name=name)

    TYPES = [Quadruped, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end
