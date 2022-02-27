"""
    RaibertHopper

    hopping robot, inspired by: "Dynamically Stable Legged Locomotion"
"""
struct RaibertHopper end

function raiberthopper(; 
    representation=:minimal, 
    timestep=0.05, 
    gravity=[0.0; 0.0; -9.81],
    control_scaling=Diagonal(ones(3)),
    seed=1, 
    contact_foot=true, 
    contact_body=true,
    vis=Visualizer(), 
    name=:robot,
    info=nothing,
    opts_step=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5), 
    opts_grad=SolverOptions(rtol=3.0e-4, btol=3.0e-4, undercut=1.5),
    T=Float64)

    mechanism = get_raiberthopper(
        timestep=timestep, 
        gravity=gravity,
        contact_foot=contact_foot,
        contact_body=contact_body)
    initialize_raiberthopper!(mechanism)

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
    control_mask = raiberthopper_control_mask()
    build_robot(mechanism, 
        vis=vis, 
        name=name)

    TYPES = [RaibertHopper, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask, control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function raiberthopper_control_mask()
    [0 0 0 1 0 0 0;
	 0 0 0 0 1 0 0;
	 0 0 0 0 0 0 1]
end

function raiberthopper_nominal_max(; 
    foot_radius=0.05)

    # initial state
    x1b1 = [0.0; 0.0; 0.5 + foot_radius]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; 0.0 + foot_radius]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function raiberthopper_offset_max(x_shift, y_shift, z_shift)
    z = raiberthopper_nominal_max()
    shift = [x_shift; y_shift; z_shift]
    z[1:3] += shift
    z[13 .+ (1:3)] += shift
    return z
end

function visualize(env::Environment{RaibertHopper}, traj::Vector{Vector{T}}; 
    name=:robot, 
    axes=false, 
    grid=true) where T

	@assert size(traj[1]) == size(env.state)

    # convert to maximal representation
    z = [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj]

    body_color = magenta
    foot_color = magenta
    env.mechanism.bodies[1].shape.color = body_color
    env.mechanism.bodies[2].shape.color = foot_color

    # build system
    build_robot(env.vis, env.mechanism, 
        name=name)

    n_leg = 100
    r_leg = 0.025
    for i = 1:n_leg
        setobject!(env.vis["leg$i"], GeometryBasics.Sphere(Point3f0(0),
            convert(Float32, r_leg)),
            MeshPhongMaterial(color=magenta))
    end

    # animate
    anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.timestep)))
    for (t, x) in enumerate(z)
        x_body = x[1:3]
        x_foot = x[13 .+ (1:3)]
        dir = x_body - x_foot
        dir_norm = dir ./ norm(dir)

        MeshCat.atframe(anim, t) do
            set_robot(env.vis, env.mechanism, x, 
                name=name)
            step = range(0.0, stop=norm(dir), length=n_leg)
            for i = 1:n_leg
                MeshCat.settransform!(env.vis["leg$i"], MeshCat.compose(MeshCat.Translation(step[i] .* dir_norm + x_foot), MeshCat.LinearMap(Rotations.RotY(0.0))))
            end
        end
    end
    MeshCat.setanimation!(env.vis, anim)

    set_camera!(env.vis, zoom=1.5)
    setvisible!(env.vis["/Axes"], axes)
    setvisible!(env.vis["/Grid"], grid)
end

function ghost(env::Environment{RaibertHopper}, traj::Vector{Vector{T}}; 
    timesteps=[t for t = 1:length(traj)], 
    axes=false, 
    grid=true, 
    line=false) where T

    # convert to maximal representation
    z = [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj]

    # color
    body_color = magenta 
    foot_color = magenta 
    env.mechanism.bodies[1].shape.color = body_color
    env.mechanism.bodies[2].shape.color = foot_color

    for t in timesteps
        # build system
        build_robot(env.vis, env.mechanism, 
            name=Symbol("robot_$t"), 
            color=(t == length(z) ? magenta : magenta_light))
        set_robot(env.vis, env.mechanism, z[t], 
            name=Symbol("robot_$t"))

        x_body = z[t][1:3]
        x_foot = z[t][13 .+ (1:3)]

        leg = GeometryBasics.Cylinder(Point3f0(x_foot...), Point3f0(x_body...), convert(Float32, 0.025))
        setobject!(env.vis["leg_$t"], leg, MeshPhongMaterial(color=(t == length(z) ? magenta : magenta_light)))
    end
    setvisible!(env.vis[:robot], false)

    # body
    if line
        line_mat = LineBasicMaterial(color=color=RGBA(51.0 / 255.0, 1.0, 1.0, 1.0), linewidth=25.0)
        points = Vector{Point{3,Float64}}()
        for (i, xt) in enumerate(z)
            k = xt[1:3]
            push!(points, Point(k[1], k[2], k[3]))
        end
        setobject!(env.vis[:body_traj], MeshCat.Line(points, line_mat))

        # foot
        line_mat = LineBasicMaterial(color=color=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0), linewidth=10.0)
        points = Vector{Point{3,Float64}}()
        for (i, xt) in enumerate(z)
            k = xt[13 .+ (1:3)]
            push!(points, Point(k[1], k[2], k[3]))
        end
        setobject!(env.vis[:foot_traj], MeshCat.Line(points, line_mat))
    end

    set_camera!(env.vis, 
        zoom=1.5)
    setvisible!(env.vis["/Axes"], axes)
    setvisible!(env.vis["/Grid"], grid)

    open(env.vis)
end
