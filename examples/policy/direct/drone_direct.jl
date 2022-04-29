# PREAMBLE

# PKG_SETUP

# ## Setup

using DirectTrajectoryOptimization
const DTO = DirectTrajectoryOptimization
using LinearAlgebra
using Plots

# ## drone
include("model_drone.jl")

nx = drone.nx
nu = drone.nu

# ## initialization
x1 = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xT = [1.0; 1.0; 0.0; 0.0; 0.0; 0.0]
θ_angle = 0.25
u_hover = hover_controls(drone, θ_angle)

# ## linear policy
# function policy(θ, x)
#     M = reshape(θ, nu, nx)
#     M * (x - xT)
# end

# nθ = nu * nx

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
T = 31

# ## model
h = 0.05

function f1(y, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    y_di = y[1:nx]
    θ = u[nu .+ (1:nθ)]
    yθ = y[nx .+ (1:nθ)]
    [
        dynamics(drone, h, y_di, x_di, u_ctrl, w);
        yθ - θ;
    ]
    # dynamics(drone, h, y_di, x_di, u_ctrl, w)
end

function ft(y, x, u, w)
    u_ctrl = u[1:nu]
    x_di = x[1:nx]
    y_di = y[1:nx]
    θ = x[nx .+ (1:nθ)]
    yθ = y[nx .+ (1:nθ)]
    [
        dynamics(drone, h, y_di, x_di, u_ctrl, w);
        yθ - θ;
    ]
    # dynamics(drone, h, y_di, x_di, u_ctrl, w)
end

dyn1 = DTO.Dynamics(f1, nx + nθ, nx, nu + nθ)
dynt = DTO.Dynamics(ft, nx + nθ, nx + nθ, nu)

dyn = [dyn1, [dynt for t = 2:T-1]...]

# ## objective
function o1(x, u, w)
    J = 0.0
    q = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    ex = x - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-2 * dot(u[1:nu] - u_hover, u[1:nu] - u_hover)
    J += 1.0e-1 * dot(u[nu .+ (1:nθ)], u[nu .+ (1:nθ)])
    return J
end

function ot(x, u, w)
    J = 0.0
    q = [10.0, 10.0, 10.0, 1.0, 1.0, 1.0]
    ex = x[1:nx] - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-2 * dot(u[1:nu] - u_hover, u[1:nu] - u_hover)
    J += 1.0e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

function ott(x, u, w)
    J = 0.0
    q = [100.0, 100.0, 100.0, 1.0, 1.0, 1.0]
    ex = x[1:nx] - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-2 * dot(u[1:nu] - u_hover, u[1:nu] - u_hover)
    J += 1.0e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

function oT(x, u, w)
    J = 0.0
    q = [10000.0, 10000.0, 10000.0, 1000.0, 1000.0, 1000.0]
    ex = x[1:nx] - xT
    J += 0.5 * transpose(ex) * Diagonal(q) * ex
    J += 1.0e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

c1 = DTO.Cost(o1, nx, nu + nθ)
ct = DTO.Cost(ot, nx + nθ, nu)
ctt = DTO.Cost(ott, nx + nθ, nu)
cT = DTO.Cost(oT, nx + nθ, 0)
obj = [c1, [ct for t = 2:16]..., [ctt for t = 17:(T - 1)]..., cT]

# ## constraints
ul = [0.0; -0.5 * π; 0.0; -0.5 * π]
uu = [10.0; 0.5 * π; 10.0; 0.5 * π]
bnd1 = DTO.Bound(nx, nu + nθ, xl=x1, xu=x1, ul=[ul; -Inf * ones(nθ)], uu=[uu; Inf * ones(nθ)])
bndt = DTO.Bound(nx + nθ, nu, ul=ul, uu=uu)
bndT = DTO.Bound(nx + nθ, 0)# xl=[xT; -Inf * ones(nθ)], xu=[xT; Inf * ones(nθ)])
bnds = [bnd1, [bndt for t = 2:T-1]..., bndT]

function policy1(x, u, w)
    θ = u[nu .+ (1:nθ)]
    u[1:nu] - policy(θ, x[1:nx], xT)
end

function policyt(x, u, w)
    θ = x[nx .+ (1:nθ)]
    u[1:nu] - policy(θ, x[1:nx], xT)
end

con_policy1 = DTO.Constraint(policy1, nx, nu + nθ)
con_policyt = DTO.Constraint(policyt, nx + nθ, nu)

cons = [con_policy1, [con_policyt for t = 2:T-1]..., DTO.Constraint()]

# ## problem
p = DTO.solver(dyn, obj, cons, bnds,
    options=DTO.Options(
        tol=1.0e-2,
        constr_viol_tol=1.0e-2))

# ## initialize
θ0 = 0.001 * randn(nθ)
u_guess = [t == 1 ? [u_hover; θ0] : u_hover for t = 2:T-1]
x_guess = [t == 1 ? x1 : [x1; randn(nθ)] for t = 1:T]

DTO.initialize_controls!(p, u_guess)
DTO.initialize_states!(p, x_guess)

# ## solve
@time DTO.solve!(p)

# ## solution
x_sol, u_sol = DTO.get_trajectory(p)
@show u_sol[1]
@show x_sol[1]
@show x_sol[T]
θ_sol = u_sol[1][nu .+ (1:nθ)]

# ## state
plot(hcat([x[1:nx] for x in x_sol]...)', label="", color=:orange, width=2.0)

# ## control
plot(hcat([u[1:nu] for u in u_sol]..., u_sol[end])', linetype = :steppost)

# ## plot xy
plot([x[1] for x in x_sol], [x[2] for x in x_sol], label="", color=:black, width=2.0)

# ## visualization
include("visuals_drone.jl")
vis = Visualizer()
open(vis)
visualize_drone!(vis, drone, [x_sol], [[u_hover, u_sol...]]; Δt=h, xT=[xT])

# ## simulate policy
x_hist = [x1]
u_hist = [u_hover]

for t = 1:(2 * T)
    push!(u_hist, policy(θ_sol, x_hist[end], [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]))
    push!(x_hist, dynamics(drone, h, x_hist[end], u_hist[end], zeros(nw)))
end

visualize_drone!(vis, drone, [x_hist], [u_hist]; Δt=h, xT=[xT])
