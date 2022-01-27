################################################################################
# Atlas
################################################################################
struct Atlas end

function atlas(; mode::Symbol=:min, dt::T=0.01, g::T=-9.81, cf=0.8,
    damper=10.0, spring=0.0, info=nothing, model_type::Symbol=:simple,
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer(), name::Symbol=:robot,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5)) where T

    mechanism = getmechanism(:atlas, timestep=dt, g=g, cf=cf, damper=damper,
        spring=spring, model_type=model_type)
    initialize!(mechanism, :atlas)

    if mode == :min
        nx = minimal_dimension(mechanism)
    elseif mode == :max
        nx = maximal_dimension(mechanism)
    end
    nu = 15
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e8 * ones(nu)), high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)
    z = get_max_state(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z
    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = [zeros(nu, 6) I(nu)]
    control_scaling = Diagonal(ones(nu))

    build_robot(vis, mechanism, name=name)

    TYPES = [Atlas, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, mode, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end
