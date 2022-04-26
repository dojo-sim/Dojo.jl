using Pkg
Pkg.activate(joinpath(module_dir(), "examples"))

using Dojo
using IterativeLQR
const iLQR = IterativeLQR
using LinearAlgebra
using Plots
# using JLD2

# ## double integrator
include("model_drone.jl")

N = 8
nx = drone.nx
nu = drone.nu

Nx = nx * N
Nu = nu * N

# ## initialization
θ_angle = 0.25
u_hover = hover_controls(drone, θ_angle)

# x1 = [
#         [-0.5; 0.125; 0.0; 0.0; 0.0; 0.0],
#         [0.25; 0.25; 0.0; 0.0; 0.0; 0.0],
#         [0.1; 0.1; 0.0; 0.0; 0.0; 0.0],
#         [-0.3; 0.05; 0.0; 0.0; 0.0; 0.0],
#         [0.6; -0.5; 0.0; 0.0; 0.0; 0.0],
#         [0.75; -0.6; 0.0; 0.0; 0.0; 0.0],
#         [-0.9; -0.2; 0.0; 0.0; 0.0; 0.0],
#         [0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
#      ]
# xT = [
#         [0.75; 0.85; 0.0; 0.0; 0.0; 0.0],
#         [0.0; 0.0; 0.0; 0.0; 0.0; 0.0],
#         [1.0; 1.0; 0.0; 0.0; 0.0; 0.0],
#         [-1.0; -1.0; 0.0; 0.0; 0.0; 0.0],
#         [-0.6; -0.5; 0.0; 0.0; 0.0; 0.0],
#         [-0.9; -1.0; 0.0; 0.0; 0.0; 0.0],
#         [-0.7; 1.0; 0.0; 0.0; 0.0; 0.0],
#         [0.85; 0.0; 0.0; 0.0; 0.0; 0.0]
#      ]

radius = 1.0
dia = sqrt(2.0) / 2.0 * radius
x1 = [
        [radius; 0.0; 0.0; 0.0; 0.0; 0.0],
        [-radius; 0.0; 0.0; 0.0; 0.0; 0.0],
        [0.0; radius; 0.0; 0.0; 0.0; 0.0],
        [0.0; -radius; 0.0; 0.0; 0.0; 0.0],
        [dia; dia; 0.0; 0.0; 0.0; 0.0],
        [-dia; -dia; 0.0; 0.0; 0.0; 0.0],
        [-dia; dia; 0.0; 0.0; 0.0; 0.0],
        [dia; -dia; 0.0; 0.0; 0.0; 0.0],
     ]
xT = [
        [-radius; 0.0; 0.0; 0.0; 0.0; 0.0],
        [radius; 0.0; 0.0; 0.0; 0.0; 0.0],
        [0.0; -radius; 0.0; 0.0; 0.0; 0.0],
        [0.0; radius; 0.0; 0.0; 0.0; 0.0],
        [-dia; -dia; 0.0; 0.0; 0.0; 0.0],
        [dia; dia; 0.0; 0.0; 0.0; 0.0],
        [dia; -dia; 0.0; 0.0; 0.0; 0.0],
        [-dia; dia; 0.0; 0.0; 0.0; 0.0],
     ]

@assert length(x1) == N
@assert length(xT) == N
X1 = vcat(x1...)
XT = vcat(xT...)

# ## (1-layer) multi-layer perceptron policy
l_input = nx
l1 = 8
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
h = 0.05

function f1(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]

    θ = u[Nu .+ (1:nθ)]

    [
        vcat([dynamics(drone, h, xi[i], ui[i], wi[i]) for i = 1:N]...)
        θ;
    ]
end

function ft(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    wi = [w[nw * (i - 1) .+ (1:nw)] for i = 1:N]

    θ = x[Nx .+ (1:nθ)]

    [
        vcat([dynamics(drone, h, xi[i], ui[i], wi[i]) for i = 1:N]...)
        θ;
    ]
end

dyn1 = iLQR.Dynamics(f1, Nx, Nu + nθ)
dynt = iLQR.Dynamics(ft, Nx + nθ, Nu)

dyn = [dyn1, [dynt for t = 2:T-1]...]

# ## objective
function o1(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]

    J = 0.0
    q = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    for i = 1:N
        ex = xi[i] - xT[i]
        eu = ui[i] - u_hover
        J += transpose(ex) * Diagonal(q) * ex  ./ N
        J += 1.0 * dot(eu, eu) ./ N
    end
    J += 1.0 * dot(θ, θ)
    return J
end

function ot(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    J = 0.0
    q = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    for i = 1:N
        ex = xi[i] - xT[i]
        eu = ui[i] - u_hover
        J += transpose(ex) * Diagonal(q) * ex  ./ N
        J += 1.0 * dot(eu, eu)
    end
    J += 1.0 * dot(θ, θ)
    return J
end

function ott(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    J = 0.0
    q = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0]
    for i = 1:N
        ex = xi[i] - xT[i]
        eu = ui[i] - u_hover
        J += transpose(ex) * Diagonal(q) * ex  ./ N
        J += 1.0 * dot(eu, eu)
    end
    J += 1.0 * dot(θ, θ)
    return J
end

function oT(x, u, w)
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = x[Nx .+ (1:nθ)]

    J = 0.0
    # q = [1000.0; 1000.0; 1500.0; 1000.0; 1000.0; 1000.0]
    # for i = 1:N
    #     # ex = xi[i] - xT[i]
    #     # J += transpose(ex) * Diagonal(q) * ex  ./ N
    # end
    J += 1.0 * dot(θ, θ)
    return J
end

c1 = iLQR.Cost(o1, Nx, Nu + nθ)
ct = iLQR.Cost(ot, Nx + nθ, Nu)
ctt = iLQR.Cost(ott, Nx + nθ, Nu)
cT = iLQR.Cost(oT, Nx + nθ, 0)
obj = [c1, [ct for t = 2:20]..., [ctt for t = 21:(T-1)]..., cT]

# ## constraints
ul = [0.0; -0.5 * π; 0.0; -0.5 * π]
uu = [100.0; 0.5 * π; 100.0; 0.5 * π]
Ul = vcat([ul for i = 1:N]...)
Uu = vcat([uu for i = 1:N]...)

function policy1(x, u, w)
    ui = [u[nu * (i - 1) .+ (1:nu)] for i = 1:N]
    xi = [x[nx * (i - 1) .+ (1:nx)] for i = 1:N]
    θ = u[Nu .+ (1:nθ)]
    [
        Ul - u[1:Nu];
        u[1:Nu] - Uu;
        vcat([ui[i] - policy(θ, xi[i], xT[i]) for i = 1:N]...)
    ]
end

function policyt(x, u, w)
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

con_policy1 = iLQR.Constraint(policy1, Nx, Nu + nθ, indices_inequality=collect(1:2Nu))
con_policyt = iLQR.Constraint(policyt, Nx + nθ, Nu, indices_inequality=collect(1:2Nu))

cons = [con_policy1, [con_policyt for t = 2:(T - 1)]..., iLQR.Constraint(goal, Nx + nθ, 0)]

# ## problem
# ## problem
opts = iLQR.Options(line_search=:armijo,
    max_iterations=250,
    max_dual_updates=10,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-1, # doesn't need to be tight for good policies
    # initial_constraint_penalty,
    # scaling_penalty,
    # max_penalty,
    # reset_cache,
    verbose=true)

# p = iLQR.problem_data(dyn, obj, cons)
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
# plt = plot();
# for i = 1:N
#     plt = plot!(hcat(X_sol[i]...)', label="", color=:orange, width=2.0)
# end
# display(plt)

# # ## control
plt = plot();
for i = 1:N
    plt = plot!(hcat(U_sol[i]..., U_sol[i][end])', linetype = :steppost)
end
display(plt)

# # ## plot xy
# plt = plot();
# for i = 1:N
#     plt = plot!([x[1] for x in X_sol[i]], [x[2] for x in X_sol[i]], label="", width=2.0)
# end
# display(plt)

# ## visualization
include("visuals_drone.jl")
vis = Visualizer()
open(vis)
visualize_drone!(vis, drone, X_sol, U_sol; Δt=h, xT=xT)

# ## simulate policy
i = 8
x_init = x1[i]
x_goal = xT[i]
x_init = [-1.0; -0.0; 0.0; 0.0; 0.0; 0.0]
x_goal = [0.25; -0.5; 0.0; 0.0; 0.0; 0.0]
x_hist = [x_init]
u_hist = [u_hover]

for t = 1:5 * T
    push!(u_hist, policy(θ_sol, x_hist[end], x_goal))
    push!(x_hist, dynamics(drone, h, x_hist[end], u_hist[end], zeros(nw)))
end

visualize_drone!(vis, drone, [x_hist], [u_hist]; Δt=h, xT=[x_goal])

# ## save policy
# @save joinpath(@__DIR__, "policy_ilqr.jld2") θ_sol
# @load joinpath(@__DIR__, "policy.jld2") θ_sol
