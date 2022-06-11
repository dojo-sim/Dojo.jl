
function newton_solver(res, x0; iterations=500, tolerance=1e-6)
    x = x0
    for i = 1:iterations
        r = res(x)
        (norm(r, Inf) < tolerance) && break
        ∇r = FiniteDiff.finite_difference_derivative(x -> res(x), x)
        Δx = - r / ∇r
        α = 1.0
        r_cand = Inf
        while norm(r_cand, Inf) >= norm(r, Inf)
            r_cand = res(x + α * Δx)
            α /= 2
        end
        x += α * Δx
        (i == iterations) && (@show "solver failure")
    end
    return x
end

function catenary_parameters(x_start, x_goal, rope_length; iterations=500, tolerance=1e-6,
        a_guess=0.0, dx_guess=0.0)

    rope_length =  max(rope_length, (1+1e-5)*norm(x_goal - x_start))
    h = x_goal[1] - x_start[1]
    v = x_goal[2] - x_start[2]

    # find a i.e. the shape of the catenary
    res(b) = 1/sqrt(sqrt(rope_length^2 - v^2)/h - 1) - 1/sqrt(2b * sinh(1/(2b)) - 1)
    b = newton_solver(res, a_guess/h, iterations=iterations, tolerance=tolerance)
    a = b * h

    # find x offset
    function caten(x; a=a)
        return a * cosh(x/a)
    end
    res = dx -> caten(x_goal[1] + dx) - caten(x_start[1] + dx) - (x_goal[2] - x_start[2])
    dx = newton_solver(res, dx_guess, iterations=iterations, tolerance=tolerance)

    # find y offset
    dy = x_goal[2] - caten(x_goal[1] + dx)

    # final function
    return a, dx, dy
end

function catenary(x; a=0.0, dx=0.0, dy=0.0)
    return a * cosh((x+dx)/a) + dy
end

function link_transform(start, goal)
    # transforms a vertical line of length 1 into a line between start and goal
    v1 = [0.0, 0.0, 1.0]
    v2 = goal[1:3] - start[1:3]
    normalize!(v2)
    ax = cross(v1, v2)
    ang = acos(v1' * v2)
    q = axis_angle_to_quaternion(ang * normalize!(ax))

    rope_length = norm(goal - start)
    scaling = LinearMap(I * Diagonal([1.0, 1.0, rope_length]))
    transform = compose(Translation(start), LinearMap(q))
    return transform, scaling
end

function build_rope(vis::Visualizer; N::Int=10, color=Colors.RGBA(0,0,0,1),
        rope_type::Symbol=:cylinder, rope_radius=0.02, name=:rope)

    if rope_type == :line
        line = MeshCat.Line([xa, xb])
        material = LineBasicMaterial(color=color, wireframeLinewidth=10, linewidth=10)
    elseif rope_type == :cylinder
        line = MeshCat.Cylinder(Point(0,0,0.0), Point(0,0,1.0), rope_radius)
        material = MeshPhongMaterial(color=color)
    end

    for i = 1:N
        setobject!(vis[name][Symbol(i)][:scaled], line, material)
    end

    return vis
end

function set_straight_rope(vis::Visualizer, x_start, x_goal; N::Int=10, name::Symbol=:rope)
    Λ = range(0,1,length=N+1)
    x = [(1-λ)*x_start + λ * x_goal for λ in Λ]
    for i = 1:N
        transform, scaling = link_transform(x[i], x[i+1])
        settransform!(vis[name][Symbol(i)][:scaled], scaling)
        settransform!(vis[name][Symbol(i)], transform)
    end
    return nothing
end

function set_loose_rope(vis::Visualizer, x_start, x_goal; N::Int=10,
        rope_length=2norm(x_goal - x_start), min_altitude=-Inf, a_guess=1.0, dx_guess=0.0, name::Symbol=:rope)
    v = x_goal - x_start
    shadow_rope_length = norm(v[1:2])
    θ = atan(v[2], v[1])
    R = rotationmatrix(RotZ(-θ))
    v̄ = R * v # rotated into the xz plane

    a, dx, dy = catenary_parameters(zeros(2), v̄[[1,3]], rope_length, a_guess=a_guess, dx_guess=dx_guess)
    Λ = shadow_rope_length * range(0,1,length=N+1)
    x = []
    for i = 1:N+1
        xi = x_start + R' * [Λ[i], 0, catenary(Λ[i], a=a, dx=dx, dy=dy)]
        xi[3] = max(xi[3], min_altitude)
        push!(x, xi)
    end
    for i = 1:N
        transform, scaling = link_transform(x[i], x[i+1])
        settransform!(vis[name][Symbol(i)][:scaled], scaling)
        settransform!(vis[name][Symbol(i)], transform)
    end
    return a, dx
end

function animate_straight_rope(vis::Visualizer, start_traj::Vector, goal_traj::Vector;
        anim::Animation=MeshCat.Animation(100), N::Int=50, name=:rope)
    M = length(start_traj)

    for i = 1:M
        atframe(anim, i) do
            xa = start_traj[i]
            xb = goal_traj[i]
            set_straight_rope(vis, xa, xb, N=N, name=name)
        end
    end
    setanimation!(vis, anim)
    return vis, anim
end

function animate_loose_rope(vis::Visualizer, start_traj::Vector, goal_traj::Vector;
        anim::Animation=MeshCat.Animation(100), rope_length=30.0, N::Int=50, min_altitude=-Inf, name=:rope)
    M = length(start_traj)

    a_guess = 1.0
    dx_guess = 0.0

    for i = 1:M
        atframe(anim, i) do
            xa = start_traj[i]
            xb = goal_traj[i]
            a_guess, dx_guess = set_loose_rope(vis, xa, xb, rope_length=rope_length,
                N=N, min_altitude=min_altitude, a_guess=a_guess, dx_guess=dx_guess, name=name)
            # set_straight_rope(vis, xa, xb, N=N, name=name)
        end
    end
    setanimation!(vis, anim)
    return vis, anim
end
