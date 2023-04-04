"""
    Block

    rigid body with fully actuation
"""
struct Block end

function block(; 
    representation=:maximal, 
    timestep=0.05, 
    gravity=-9.81,
    friction_coefficient=0.8, 
    side=0.5, 
    contact=true, 
    contact_type=:nonlinear,
    seed=1, 
    vis=Visualizer(), 
    info=nothing, 
    name=:robot,
    control_scaling=Diagonal(ones(3)),
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5), 
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_block(;
        timestep, 
        gravity, 
        friction_coefficient, 
        side, 
        contact, 
        contact_type)

    initialize_block!(mechanism)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu = 3
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
    control_mask = [I(nu) zeros(nu, nu)]

    build_robot(mechanism, 
        vis=vis, 
        name=name)

    TYPES = [Block, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function visualize_block_force!(vis, anim, z_vis, u_vis; 
    shift = [0; 0.25; 0], 
    color=orange,
    )

    u_max = maximum([norm(u) for u in u_vis])
    force_vis = Dojo.ArrowVisualizer(vis[:force])
    setobject!(force_vis, Dojo.MeshPhongMaterial(color=color))
  
    for t = 1:length(z_vis)
        z = z_vis[t]
        u = (t == length(z_vis) ? 0.5 * u_vis[end] ./ u_max : 0.5 * u_vis[t] ./ u_max)
        MeshCat.atframe(anim, t) do
            settransform!(force_vis,
                Dojo.Point(z[1] - u[1] - shift[1], z[2] - u[2] - shift[2], z[3] - u[3] - shift[3]),
                Dojo.Vec(u[1], u[2], u[3]),
                shaft_radius=0.05,
                max_head_radius=0.1)
        end 
    end
    MeshCat.setanimation!(vis, anim)

    return vis, anim
end