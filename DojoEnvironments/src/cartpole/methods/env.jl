"""
    Cartpole <: Environment

    classic system with one translational (cart) and one rotational (pole) degree of freedom
    https://underactuated.mit.edu/acrobot.html#cart_pole
"""
struct Cartpole end

function cartpole(; 
    representation=:minimal, 
    timestep=0.05, 
    gravity=[0.0; 0.0; -9.81],
    spring=0.0, 
    damper=0.0,
    seed=1, 
    vis=Visualizer(), 
    info=nothing, 
    name=:robot,
    control_scaling=Diagonal(ones(1)),
    opts_step=SolverOptions(), 
    opts_grad=SolverOptions(),
    T=Float64)

    mechanism = get_cartpole(;
        timestep, 
        gravity,
        spring,
        damper)
    initialize_cartpole!(mechanism)

    if representation == :minimal
        nx = minimal_dimension(mechanism)
    elseif representation == :maximal
        nx = maximal_dimension(mechanism)
    end
    nu = 1
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
    control_mask = [1.0 0.0]

    build_robot(mechanism, vis=vis, name=name)

    TYPES = [Cartpole, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, representation, aspace, ospace,
        x, fx, fu,
        u_prev, control_mask' * control_scaling,
        nx, nu, no,
        info,
        [rng], vis,
        opts_step, opts_grad)

    return env
end

function cartpole_nominal_max(;
    pendulum_length=1.0)

    # initial state
    x1b1 = [0.0; 0.0; 0.0]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; -0.5 * pendulum_length]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function cartpole_goal_max(;
    pendulum_length=1.0)

    # target state
    xTb1 = [0.0; 0.0; 0.0]
    vTb1 = [0.0; 0.0; 0.0]
    qTb1 = [1.0; 0.0; 0.0; 0.0]
    ωTb1 = [0.0; 0.0; 0.0]
    zTb1 = [xTb1; vTb1; qTb1; ωTb1]

    xTb2 = [0.0; 0.0; 0.5 * pendulum_length]
    vTb2 = [0.0; 0.0; 0.0]
    qTb2 = [0.0; 1.0; 0.0; 0.0]
    ωTb2 = [0.0; 0.0; 0.0]
    zTb2 = [xTb2; vTb2; qTb2; ωTb2]

    zT = [zTb1; zTb2]
end

function visualize(env::Environment{Cartpole}, traj::Vector{Vector{T}}; 
    axes=false, 
    grid=false, 
    ee_traj=false) where T

	@assert size(traj[1]) == size(env.state)
    storage = generate_storage(env.mechanism, [env.representation == :minimal ? minimal_to_maximal(env.mechanism, x) : x for x in traj])
    visualize(env.mechanism, storage, 
        vis=env.vis)
	set_camera!(env.vis, zoom=1.0, cam_pos=[2,0,0])
    setvisible!(env.vis["/Axes"], axes)
    setvisible!(env.vis["/Grid"], grid)
    # slider
    l2 = GeometryBasics.Cylinder(Point3f0(0.0, -5.0, 0.0),
		Point3f0(0.0, 5.0, 0.0),
		convert(Float32, 0.0125))
	setobject!(env.vis["slider"], l2, MeshPhongMaterial(color=Colors.RGBA(0.0, 0.0, 0.0, 1.0)))
end

function ghost(env::Environment{Cartpole}, z_sol; 
    timesteps=[t for t = 1:length(z_sol)], 
    line=false)

    for t in timesteps
        name = Symbol("robot_$t")
        build_robot(env.vis, env.mechanism, 
            name=name, 
            color=(t == length(z_sol) ? cyan : cyan_light))
        z = z_sol[t]
        if t == length(z_sol)
            z[1] -= 0.05
            z[14] -= 0.05
        end
        set_robot(env.vis, env.mechanism, z, 
            name=name)
    end

	set_camera!(env.vis, 
        zoom=1.0, 
        cam_pos=[2,0,0])
    setvisible!(env.vis[:robot], false)
    l2 = GeometryBasics.Cylinder(Point3f0(0.0, -5.0, 0.0),
            Point3f0(0.0, 5.0, 0.0),
            convert(Float32, 0.0125))
    setobject!(env.vis["slider"], l2, MeshPhongMaterial(color=Colors.RGBA(0.0, 0.0, 0.0,1.0)))
    setvisible!(env.vis["/Axes"], false)
    setvisible!(env.vis["/Grid"], false)

    if line
        line_mat = LineBasicMaterial(
            color=color=RGBA(51.0 / 255.0, 1.0, 1.0, 1.0), 
            linewidth=25.0)
        points = Vector{Point{3,Float64}}()
        for (i, xt) in enumerate(z_sol)
            z_pendulum = xt[13 .+ (1:13)]
            x = z_pendulum[1:3]
            q = Quaternion(z_pendulum[6 .+ (1:4)]...)
            k = x + vector_rotate([0.0; 0.0; -0.5], q)

            push!(points, Point(k[1], k[2], k[3]))

            setobject!(env.vis["ee_vertex_$i"], GeometryBasics.Sphere(Point3f0(0),
                convert(Float32, 0.001)),
                MeshPhongMaterial(color=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0)))
                settransform!(env.vis["ee_vertex_$i"], MeshCat.Translation(points[i]))
        end
        setobject!(env.vis[:ee_traj], MeshCat.Line(points, line_mat))
    end

    open(env.vis)
end
