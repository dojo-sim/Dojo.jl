################################################################################
# RaibertHopper
################################################################################
struct RaibertHopper end

function raiberthopper(; mode::Symbol=:min, dt::T=0.05, g::T=-9.81,
    control_scaling=Diagonal(ones(3)),
    s::Int=1, contact::Bool=true, vis::Visualizer=Visualizer(), name::Symbol=:robot,
    info=nothing,
    opts_step=InteriorPointOptions(), opts_grad=InteriorPointOptions()) where T

    mechanism = getraiberthopper(Δt=dt, g=g)
    initializeraiberthopper!(mechanism)

    if mode == :min
        nx = minCoordDim(mechanism)
    elseif mode == :max
        nx = maxCoordDim(mechanism)
    end
    nu = 3
    no = nx

    aspace = BoxSpace(nu, low=(-1.0e8 * ones(nu)), high=(1.0e8 * ones(nu)))
    ospace = BoxSpace(no, low=(-Inf * ones(no)), high=(Inf * ones(no)))

    rng = MersenneTwister(s)

    z = getMaxState(mechanism)
    x = mode == :min ? max2min(mechanism, z) : z

    fx = zeros(nx, nx)
    fu = zeros(nx, nu)

    u_prev = zeros(nu)
    control_mask = raiberthopper_control_mask()
    build_robot(vis, mechanism, name=name)

    TYPES = [RaibertHopper, T, typeof(mechanism), typeof(aspace), typeof(ospace), typeof(info)]
    env = Environment{TYPES...}(mechanism, mode, aspace, ospace,
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

function raiberthopper_nominal_max(; foot_radius=0.05)
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

function visualize(env::Environment{RaibertHopper}, traj::Vector{Vector{T}}; name=:robot, axes=false, grid=true) where T
    # convert to maximal representation
    z = [env.mode == :min ? min2max(env.mechanism, x) : x for x in traj]

    body_color = orange#RGBA(0.0, 0.0, 0.0, 1.0)
    foot_color = orange#RGBA(0.0, 0.0, 0.0, 1.0) #RGBA(1.0, 165.0 / 255.0, 0.0, 1.0) 
    env.mechanism.bodies[3].shape.color = body_color 
    env.mechanism.bodies[4].shape.color = foot_color 

    # build system
    build_robot(env.vis, env.mechanism, name=name)

    n_leg = 100
    r_leg = 0.025
    for i = 1:n_leg
        setobject!(env.vis["leg$i"], GeometryBasics.Sphere(Point3f0(0),
            convert(Float32, r_leg)),
            MeshPhongMaterial(color=orange))
    end

    # animate
    anim = MeshCat.Animation(convert(Int, floor(1.0 / env.mechanism.Δt)))
    for (t, x) in enumerate(z) 
        x_body = x[1:3] 
        x_foot = x[13 .+ (1:3)] 
        dir = x_body - x_foot 
        dir_norm = dir ./ norm(dir) 

        MeshCat.atframe(anim, t) do
            set_robot(env.vis, env.mechanism, x, name=name)
            step = range(0.0, stop=norm(dir), length=n_leg) 
            for i = 1:n_leg 
                MeshCat.settransform!(env.vis["leg$i"], MeshCat.compose(MeshCat.Translation(step[i] .* dir_norm + x_foot), MeshCat.LinearMap(Rotations.RotY(0.0))))
            end
        end
    end
    MeshCat.setanimation!(env.vis, anim)

    # MeshCat.settransform!(env.vis["/Cameras/default"],
    #     MeshCat.compose(MeshCat.Translation(0.0, 0.0, -1.0), MeshCat.LinearMap(Rotations.RotZ(-1.0 * pi))))
    setprop!(env.vis["/Cameras/default/rotated/<object>"], "zoom", 1.5)
    setvisible!(env.vis["/Axes"], axes)
    setvisible!(env.vis["/Grid"], grid)
end

function ghost(env::Environment{RaibertHopper}, traj::Vector{Vector{T}}; timesteps=[t for t = 1:length(traj)], axes=false, grid=true, line=false) where T
    # convert to maximal representation
    z = [env.mode == :min ? min2max(env.mechanism, x) : x for x in traj]

    # color 
    body_color = orange #RGBA(0.0, 0.0, 0.0, 1.0)
    foot_color = orange #RGBA(0.0, 0.0, 0.0, 1.0)#RGBA(1.0, 165.0 / 255.0, 0.0, 1.0) 
    env.mechanism.bodies[3].shape.color = body_color 
    env.mechanism.bodies[4].shape.color = foot_color 

    for t in timesteps
        # build system
        build_robot(env.vis, env.mechanism, name=Symbol("robot_$t"), color=(t == length(z) ? orange : orange_light))
        set_robot(env.vis, env.mechanism, z[t], name=Symbol("robot_$t"))

        x_body = z[t][1:3] 
        x_foot = z[t][13 .+ (1:3)] 

        leg = GeometryBasics.Cylinder(Point3f0(x_foot...), Point3f0(x_body...), convert(Float32, 0.025))
        setobject!(env.vis["leg_$t"], leg, MeshPhongMaterial(color=(t == length(z) ? orange : orange_light)))
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

    # MeshCat.settransform!(env.vis["/Cameras/default"],
    #     MeshCat.compose(MeshCat.Translation(0.0, 0.0, -1.0), MeshCat.LinearMap(Rotations.RotZ(-1.0 * pi))))
    setprop!(env.vis["/Cameras/default/rotated/<object>"], "zoom", 1.5)
    setvisible!(env.vis["/Axes"], axes)
    setvisible!(env.vis["/Grid"], grid)

    open(env.vis)
end