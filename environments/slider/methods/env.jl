################################################################################
# Slider
################################################################################
struct Slider end

function slider(;
    representation=:minimal,
    timestep=0.05,
    gravity=[0.0; 0.0; 0.0],
    spring=0.0,
    damper=0.0,
    info=nothing,
    seed=1,
    vis=Visualizer(),
    name=:robot,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_mechanism(:slider,
        timestep=timestep,
        gravity=gravity,
        spring=spring,
        damper=damper)

    initialize!(mechanism, :slider)

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

    TYPES = [Slider, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end
