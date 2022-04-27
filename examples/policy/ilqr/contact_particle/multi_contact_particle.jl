using Pkg
Pkg.activate(joinpath(module_dir(), "examples"))

using Dojo
using IterativeLQR
const iLQR = IterativeLQR
using LinearAlgebra
using Plots
using Symbolics
using BenchmarkTools


# ## contact particle
vis = Visualizer()
open(vis)


include("model_contact_particle.jl")

timestep = 0.05
gravity = -9.81
rtol = 1e-7
btol = 1e-5
model = ContactParticle19(gravity=gravity, timestep=timestep)
add_symbolics!(model)
env = particle(timestep=timestep)

N = 8
nx = model.nx
nu = model.nu
nw = 0
Nx = nx * N
Nu = nu * N

# ## initialization
radius = 1.0
dia = sqrt(2.0) / 2.0 * radius
x1 = [
        [+radius; 0.0; 0.25; 0.0; 0.0; 0.0],
        [-radius; 0.0; 0.25; 0.0; 0.0; 0.0],
        [0.0; +radius; 0.25; 0.0; 0.0; 0.0],
        [0.0; -radius; 0.25; 0.0; 0.0; 0.0],
        [+dia; +dia; 0.25; 0.0; 0.0; 0.0],
        [-dia; -dia; 0.25; 0.0; 0.0; 0.0],
        [-dia; +dia; 0.25; 0.0; 0.0; 0.0],
        [+dia; -dia; 0.25; 0.0; 0.0; 0.0],
     ]
xT = [
        [-radius; 0.0; 0.25; 0.0; 0.0; 0.0],
        [+radius; 0.0; 0.25; 0.0; 0.0; 0.0],
        [0.0; -radius; 0.25; 0.0; 0.0; 0.0],
        [0.0; +radius; 0.25; 0.0; 0.0; 0.0],
        [-dia; -dia; 0.25; 0.0; 0.0; 0.0],
        [+dia; +dia; 0.25; 0.0; 0.0; 0.0],
        [+dia; -dia; 0.25; 0.0; 0.0; 0.0],
        [-dia; +dia; 0.25; 0.0; 0.0; 0.0],
     ]

u_hover = [0.0, 0.0, 0.0]

@assert length(x1) == N
@assert length(xT) == N
X1 = vcat(x1...)
XT = vcat(xT...)


# ## (1-layer) multi-layer perceptron policy
l_input = nx
l1 = 6
l2 = nu
nθ = l1 * l_input + l2 * l1

function policy(θ, x, goal)
    shift = 0
    # input
    input = x - goal

    # layer 1
    W1 = reshape(θ[shift .+ (1:(l1 * l_input))], l1, l_input)
    z1 = W1 * input
    o1 = tanh.(z1)
    shift += l1 * l_input

    # layer 2
    W2 = reshape(θ[shift .+ (1:(l2 * l1))], l2, l1)
    z2 = W2 * o1

    o2 = z2
    return o2
end

# ## horizon
T = 41

# ## model
h = timestep

function f1(y, x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]

    for i = 1:N
        dynamics(model, view(y, (i-1)*nx .+ (1:nx)), xi[i], ui[i])
    end
    y[Nx .+ (1:nθ)] .= θ
end

function f1x(dx, x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]

    dx .= 0.0
    for i = 1:N
        dynamics_jacobian_state(model, view(dx, (i-1)*nx .+ (1:nx), (i-1)*nx .+ (1:nx)), xi[i], ui[i])
    end
end

function f1u(du, x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]

    du .= 0.0
    for i = 1:N
        dynamics_jacobian_input(model, view(du, (i-1)*nx .+ (1:nx), (i-1)*nu .+ (1:nu)), xi[i], ui[i])
    end
    du[Nx .+ (1:nθ), Nu .+ (1:nθ)] .= I(nθ)
end


function ft(y, x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    for i = 1:N
        dynamics(model, view(y, (i-1)*nx .+ (1:nx)), xi[i], ui[i])
    end
    y[Nx .+ (1:nθ)] .= θ
end

function ftx(dx, x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    dx .= 0.0
    for i = 1:N
        dynamics_jacobian_state(model, view(dx, (i-1)*nx .+ (1:nx), (i-1)*nx .+ (1:nx)), xi[i], ui[i])
    end
    dx[Nx .+ (1:nθ), Nx .+ (1:nθ)] .= I(nθ)
end

function ftu(du, x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    for i = 1:N
        dynamics_jacobian_input(model, view(du, (i-1)*nx .+ (1:nx), (i-1)*nu .+ (1:nu)), xi[i], ui[i])
    end
end



# user-provided dynamics and gradients
dyn1 = iLQR.Dynamics(f1, f1x, f1u, Nx + nθ, Nx, Nu + nθ)
dynt = iLQR.Dynamics(ft, ftx, ftu, Nx + nθ, Nx + nθ, Nu)

dyn = [dyn1, [dynt for t = 2:T-1]...]

# ## objective
function o1(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]

    J = 0.0
    q = [1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1]
    r = [1.0e-2, 1.0e-2, 1.0e-0]
    for i = 1:N
        ex = xi[i] - xT[i]
        eu = ui[i] - u_hover
        J += 0.5 * transpose(ex) * Diagonal(q) * ex  ./ N
        J += 0.5 * transpose(eu) * Diagonal(r) * eu  ./ N
    end
    J += 1e-1 * dot(θ, θ)
    return J
end

function ot(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    J = 0.0
    q = [1.0, 1.0, 1.0, 1.0e-1, 1.0e-1, 1.0e-1]
    r = [1.0e-2, 1.0e-2, 1.0e-0]
    for i = 1:N
        ex = xi[i] - xT[i]
        eu = ui[i] - u_hover
        J += 0.5 * transpose(ex) * Diagonal(q) * ex  ./ N
        J += 0.5 * transpose(eu) * Diagonal(r) * eu  ./ N
    end
    J += 1e-1 * dot(θ, θ)
    return J
end

function ott(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    J = 0.0
    q = [10.0, 10.0, 10.0, 1.0e-1, 1.0e-1, 1.0e-1]
    r = [1.0e-2, 1.0e-2, 1.0e-0]
    for i = 1:N
        ex = xi[i] - xT[i]
        eu = ui[i] - u_hover
        J += transpose(ex) * Diagonal(q) * ex  ./ N
        J += 0.5 * transpose(eu) * Diagonal(r) * eu  ./ N
    end
    J += 1e-1 * dot(θ, θ)
    return J
end

function oT(x, u, w)
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    J = 0.0
    J += 1e-1 * dot(θ, θ)
    return J
end

c1 = iLQR.Cost(o1, Nx, Nu + nθ)
ct = iLQR.Cost(ot, Nx + nθ, Nu)
ctt = iLQR.Cost(ott, Nx + nθ, Nu)
cT = iLQR.Cost(oT, Nx + nθ, 0)
obj = [c1, [ct for t = 2:16]..., [ctt for t = 17:(T - 1)]..., cT]

# ## constraints
Ul = -10.0 * ones(Nu)
Uu = +10.0 * ones(Nu)

function con1(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]
    [
        Ul - u[1:Nu];
        u[1:Nu] - Uu;
        vcat([ui[i] - policy(θ, xi[i], xT[i]) for i = 1:N]...)
    ]
end

function cont(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]
    [
        Ul - u[1:Nu];
        u[1:Nu] - Uu;
        vcat([ui[i] - policy(θ, xi[i], xT[i]) for i = 1:N]...)
    ]
end

function goal(x, u, w)
    x[1:Nx] - XT[1:Nx]
end
con_policy1 = iLQR.Constraint(con1, Nx, Nu + nθ, indices_inequality=collect(1:2Nu))
con_policyt = iLQR.Constraint(cont, Nx + nθ, Nu, indices_inequality=collect(1:2Nu))

cons = [con_policy1, [con_policyt for t = 2:T-1]..., iLQR.Constraint(goal, Nx + nθ, 0)]

# ## problem
opts = iLQR.Options(line_search=:armijo,
    max_iterations=250,
    max_dual_updates=10,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-3,
    verbose=true)

p = iLQR.Solver(dyn, obj, cons, options=opts)

# ## initialize
θ0 = 1.0 * randn(nθ)
U_hover = vcat([u_hover for i = 1:N]...)
u_guess = [[U_hover; θ0], [U_hover for t = 2:T-1]...]
x_guess = iLQR.rollout(dyn, X1, u_guess)

iLQR.initialize_controls!(p, u_guess)
iLQR.initialize_states!(p, x_guess)

# ## solve
@time iLQR.solve!(p)

# ## solution
x_sol, u_sol = iLQR.get_trajectory(p)
X_sol = [[x[(i - 1) * nx .+ (1:nx)] for x in x_sol] for i = 1:N]
U_sol = [[u_hover, [u[(i - 1) * nu .+ (1:nu)] for u in u_sol]...] for i = 1:N]
@show θ_sol = u_sol[1][Nu .+ (1:nθ)]
@show x_sol[2][Nx .+ (1:nθ)]
θ_sol - x_sol[2][Nx .+ (1:nθ)]



# # ## state
plt = plot();
for i = 1:N
    plt = plot!(hcat(X_sol[i]...)', label="", color=:orange, width=2.0)
end
display(plt)

# # ## control
plt = plot();
for i = 1:N
    plt = plot!(hcat(U_sol[i]..., U_sol[i][end])', linetype = :steppost)
end
display(plt)

# ## plot xy
plt = plot();
for i = 1:N
    plt = plot!([x[1] for x in X_sol[i]], [x[2] for x in X_sol[i]], label="", width=2.0)
end
display(plt)


# ## visualization
Z_sol = [[minimal_to_maximal(env.mechanism, x) for x in X_sol[i]] for i=1:N]
storage_sol = [generate_storage(env.mechanism, Z_sol[i]) for i=1:N]
for i = 1:N
    vis, anim = visualize(env.mechanism, storage_sol[i], vis=vis, animation=anim, name=Symbol(:robot_, i))
end


# ## simulate policy
i = 8
x_init = x1[i]
x_goal = xT[i]
x_init = [-3.0; -1.0; 0.25; 0.0; 0.0; 0.0]
x_goal = [0.0; -0.0; 0.25; 0.0; 0.0; 0.0]
x_hist = [x_init]
u_hist = [u_hover]

for t = 1:3 * T
    push!(u_hist, policy(θ_sol, x_hist[end], x_goal))
    y = zeros(nx)
    dynamics(model, y, x_hist[end], u_hist[end])
    push!(x_hist, y)
end

z_hist = [minimal_to_maximal(env.mechanism, x) for x in x_hist]
storage_hist = generate_storage(env.mechanism, z_hist)
visualize(env.mechanism, storage_hist, vis=vis)

# convert_frames_to_video_and_gif("multi_particle_policy_generalization")
# convert_frames_to_video_and_gif("multi_particle_open_loop")
