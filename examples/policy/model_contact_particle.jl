struct Particle{T}
    nx::Int
    nu::Int
    nw::Int
    mass::T
    gravity::T
    side_length::T
    timestep::T
end

function residual(model::Particle, x, u, w)
    v15 = x[4:6]
    p2 = x[1:3]
    p1 = p2 - model.timestep * v15

    dyn = model.mass * ()
    res = [dyn]
    return res
end

nx = 6
nu = 3
nw = 3
model = Particle(nx, nu, nw, 1.0, -9.81, 0.5)
x = zeros(nx)
u = zeros(nu)
w = zeros(nw)
residual(model, x, u, w)

# continuous time dynamics
function dynamics(model::Particle, x, u, w)
    mass_matrix = Diagonal([model.mass, model.mass, model.mass])
    dynamics_bias = [0.0, 0.0, -model.mass * model.gravity]

    ## control mapping
    # u = (t1, θ1, t2, θ2)
    f1 = rotation_matrix(u[2]) * [0.0; u[1]]
    f2 = rotation_matrix(u[4]) * [0.0; u[3]]

    r1 = [0.5 * model.side_length; 0.0]
    m1 = r1[1] * f1[2] # - r1[2] * f1[1]

    r2 = [-0.5 * model.side_length; 0.0]
    m2 = r2[1] * f2[2] # - r2[2] * f2[1]

    rot_body = rotation_matrix(x[3])

    input = [rot_body * f1; m1] + [rot_body * f2; m2]

    # dynamics
    qdd = mass_matrix \ (dynamics_bias + input)

    return [x[4:6]; qdd]
end

# discrete-time dynamics (implicit midpoint)
function dynamics(model::Particle, h, y, x, u, w)
    y - (x + h * dynamics(model, 0.5 * (x + y), u, w))
end

# discrete-time dynamics (explicit midpoint)
function dynamics(model::Particle, h, x, u, w)
    x + h * dynamics(model, x + 0.5 * h * dynamics(model, x, u, w), u, w)
end

function hover_controls(drone::Particle, θ)
    [0.5 * drone.mass * drone.gravity / cos(θ); θ; 0.5 * drone.mass * drone.gravity / cos(θ); -θ]
end

# model
nx = 6
nu = 4
nw = 0
mass = 1.0
side_length = 0.5
inertia = 1.0 / 12.0 * mass * side_length^2
gravity = 9.81
drone = Particle(nx, nu, nw, mass, inertia, gravity, side_length)

# u_hover = [0.5 * drone.mass * drone.gravity; 0.0; 0.5 * drone.mass * drone.gravity; 0.0]

# θ = 0.35
# u_hover = hover_controls(drone, θ)

# T = 10
# h = 0.1
# x_hist = [[0.1, 0.0, -0.3, 0.0, 0.0, 0.0]]
# u_hist = Vector{eltype(drone.mass)}[]

# for t = 1:T
#     push!(u_hist, u_hover)
#     push!(x_hist, dynamics(drone, h, x_hist[end], u_hist[end], zeros(drone.nw)))
# end

# vis = Visualizer()
# render(vis)

# visualize!(vis, drone, [x_hist], [u_hist], Δt=h)
