################################################################################
# Nerf & Sphere
################################################################################
struct NerfSphere end

function nerf_sphere(;
    representation=:minimal,
    nerf::Symbol=:bunny,
    collider_options=ColliderOptions(),
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.5,
    mass=10.0,
    info=nothing,
    seed=1,
    contact=true,
    vis=Visualizer(),
    name=:robot,
    infeasible_control=false,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_mechanism(:nerf_sphere,
        nerf=nerf,
        collider_options=collider_options,
        timestep=timestep,
        gravity=gravity,
        friction_coefficient=friction_coefficient,
        radius=radius,
        mass=mass)

    initialize!(mechanism, :nerf_sphere)

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

    TYPES = [NerfSphere, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)
end
