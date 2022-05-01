using Pkg
Pkg.activate(joinpath(Dojo.module_dir(), "examples"))

using Dojo
using IterativeLQR
using RoboDojo
using Plots
using Symbolics
using BenchmarkTools
using LinearAlgebra
using FiniteDiff
using StaticArrays

const iLQR = IterativeLQR
const RD = RoboDojo

include("../methods.jl")

vis = Visualizer()
open(vis)

include("gait_design.jl")

include("../../robodojo/centroidal_quadruped/model.jl")
include("../../robodojo/centroidal_quadruped/visuals.jl")
include("../../robodojo/centroidal_quadruped/simulator.jl")
include("../../robodojo/dynamics.jl")

RoboDojo.RESIDUAL_EXPR
force_codegen = true
# force_codegen = false
robot = centroidal_quadruped
include("../../robodojo/codegen.jl")
RoboDojo.RESIDUAL_EXPR


################################################################################
# Recover data
################################################################################
file = JLD2.jldopen(joinpath(@__DIR__, "centroidal_quadruped_sol.jld2"))
x_ref = file["x_sol"]
x_sol = file["x_sol"]
u_sol = file["u_sol"]
K_sol = file["K_sol"]


################################################################################
# Policy  & Derivatives
################################################################################

# ## (2-layer) multi-layer perceptron policy
l_input = nx
l1 = 12
nθ = l1 * l_input

function feature(x)
    xbody = x[1:3]
    xorient = x[4:6]
    xfoot1 = x[7:9]
    xfoot2 = x[10:12]
    xfoot3 = x[13:15]
    xfoot4 = x[16:18]
    xf = [xbody[2:3]; xorient; xfoot1 - xbody; xfoot2 - xbody; xfoot3 - xbody; xfoot4 - xbody; x[nq .+ (1:nq)]]
    return xf
end

function policy(θ, x, w)
    shift = 0
    # layer 0
    input = [1; feature(x) - feature(w[1:nx])] # policy independent of the x position

    # layer 1
    W1 = reshape(θ[shift .+ (1:(l1 * l_input))], l1, l_input)
    z1 = W1 * input
    return z1
end

@variables x_[1:nx]
@variables θ_[1:nθ]
@variables w_[1:nx]

pi = policy(θ_, x_, w_)
∂pi∂x = Symbolics.jacobian(pi, x_)
∂pi∂xv = reshape(∂pi∂x, (nu-nu_infeasible)*nx)
∂pi∂θ = Symbolics.jacobian_sparsity(pi, θ_)
∂pi∂xvθ = Symbolics.jacobian_sparsity(∂pi∂xv , θ_)

f_pi = eval(build_function(pi, θ_, x_, w_, checkbounds=true)[2])
f_∂pi∂x = eval(build_function(∂pi∂x, θ_, x_, w_, checkbounds=true)[2])
f_∂pi∂θ = eval(build_function(∂pi∂θ, θ_, x_, w_, checkbounds=true)[2])
f_∂pi∂xv = eval(build_function(∂pi∂xv, θ_, x_, w_, checkbounds=true)[2])
f_∂pi∂xvθ = eval(build_function(∂pi∂xvθ, θ_, x_, w_, checkbounds=true)[2])
∂pi∂x0 = similar(∂pi∂x, Float64) # zeros(nu-nu_infeasible, nx)
∂pi∂θ0 = similar(∂pi∂θ, Float64) # spzeros(nu-nu_infeasible, nθ)
∂pi∂xv0 = similar(∂pi∂xv, Float64)
∂pi∂xvθ0 = similar(∂pi∂xvθ, Float64) # spzeros((nu-nu_infeasible) * nx, nθ)
pi0 = zeros(nu-nu_infeasible)
θ0 = zeros(nθ)
x0 = zeros(nx)
w0 = zeros(nx)

# @benchmark $f_∂pi∂xvθ($∂pi∂xvθ0, $θ0, $x0, $w0)
# @benchmark $f_∂pi∂xv($∂pi∂xv0, $θ0, $x0, $w0)
# @benchmark $f_∂pi∂θ($∂pi∂θ0, $θ0, $x0, $w0)
# @benchmark $f_∂pi∂x($∂pi∂x0, $θ0, $x0, $w0)
# @benchmark $f_pi($pi0, $θ0, $x0, $w0)

mutable struct PolicyCache113{T}
    c::T
    g::Vector{T}
    H::SparseMatrixCSC{T, Int}

    pi::Vector{T}
    ∂pi∂θ::SparseMatrixCSC{T, Int}
    ∂pi∂xv::Vector{T}
    ∂pi∂xvθ::SparseMatrixCSC{T, Int}

    Qu::Diagonal{T, Vector{T}}
    QKv::Diagonal{T, Vector{T}}
end

function PolicyCache113(nu, nu_infeasible, nx, nθ, ∂pi∂θ_sparsity, ∂pi∂xvθ_sparsity)
    c = 0.0
    g = zeros(nθ)
    H = spzeros(nθ, nθ)

    pi = zeros(nu-nu_infeasible)
    ∂pi∂θ = similar(∂pi∂θ_sparsity, Float64)
    ∂pi∂xv = zeros((nu-nu_infeasible)*nx)
    ∂pi∂xvθ = similar(∂pi∂xvθ_sparsity, Float64)

    Qu = 1e0 .* Diagonal(ones(nu-nu_infeasible))
    QKv = 1e-2 .* Diagonal(ones((nu-nu_infeasible)*nx))

    return PolicyCache113(
        c, g, H,
        pi,
        ∂pi∂θ,
        ∂pi∂xv,
        ∂pi∂xvθ,
        Qu, QKv,
        )
end

function evaluation_and_gradient!(policy_cache::PolicyCache113, data::LQRData113, θ)
    pi = policy_cache.pi
    ∂pi∂θ = policy_cache.∂pi∂θ
    ∂pi∂xv = policy_cache.∂pi∂xv
    ∂pi∂xvθ = policy_cache.∂pi∂xvθ
    Qu = policy_cache.Qu
    QKv = policy_cache.QKv
    T = length(x_ref)

    policy_cache.c = 0.0
    policy_cache.g .= 0.0
    policy_cache.H .= 0.0

    for i = 1:T-1
        xi = data.x_sol[i]
        ui = data.u_sol[i][nu_infeasible+1:nu]
        Kvi = reshape(data.K_sol[i][nu_infeasible+1:nu, :], (nu-nu_infeasible)*nx)
        wi = data.x_ref[i]

        f_pi(policy_cache.pi, θ, xi, wi)
        f_∂pi∂θ(policy_cache.∂pi∂θ, θ, xi, wi)
        f_∂pi∂xv(policy_cache.∂pi∂xv, θ, xi, wi)
        f_∂pi∂xvθ(policy_cache.∂pi∂xvθ, θ, xi, wi)

        policy_cache.c += 0.5 * (pi - ui)' * Qu *(pi - ui)
        policy_cache.c += 0.5 * (∂pi∂xv - Kvi)' * QKv *(∂pi∂xv - Kvi)

        policy_cache.g += ∂pi∂θ' * Qu * (pi - ui)
        policy_cache.g += ∂pi∂xvθ' * QKv * (∂pi∂xv - Kvi)

        policy_cache.H += ∂pi∂θ' * Qu * ∂pi∂θ
        policy_cache.H += ∂pi∂xvθ' * QKv * ∂pi∂xvθ
    end
    return nothing
end

function evaluation!(policy_cache::PolicyCache113, data::LQRData113, θ)
    T = length(x_ref)
    pi = policy_cache.pi
    ∂pi∂xv = policy_cache.∂pi∂xv
    Qu = policy_cache.Qu
    QKv = policy_cache.QKv

    policy_cache.c = 0.0

    for i = 1:T-1
        xi = data.x_sol[i]
        ui = data.u_sol[i][nu_infeasible+1:nu]
        Kvi = reshape(data.K_sol[i][nu_infeasible+1:nu, :], (nu-nu_infeasible)*nx)
        wi = data.x_ref[i]

        f_pi(policy_cache.pi, θ, xi, wi)
        f_∂pi∂xv(policy_cache.∂pi∂xv, θ, xi, wi)

        policy_cache.c += 0.5 * (pi - ui)' * Qu *(pi - ui)
        policy_cache.c += 0.5 * (∂pi∂xv - Kvi)' * QKv *(∂pi∂xv - Kvi)

    end
    return nothing
end

mutable struct LQRData113{T}
    x_ref::Vector{Vector{T}}
    x_sol::Vector{Vector{T}}
    u_sol::Vector{Vector{T}}
    K_sol::Vector{Matrix{T}}
end

data = LQRData113(x_ref, x_sol, u_sol, K_sol)
policy_cache = PolicyCache113(nu, nu_infeasible, nx, nθ, ∂pi∂θ, ∂pi∂xvθ)

# @benchmark evaluation_and_gradient!(policy_cache, data, θ0)
# @benchmark evaluation!(policy_cache, data, θ0)


function newton_solve(policy_cache::PolicyCache113, data::LQRData113)
    θ = zeros(nθ)

    for i = 1:10
        (c <= 1e-3) && break
        evaluation_and_gradient!(policy_cache, θ, data)
        Δ = - (policy_cache.H + I) \ policy_cache.g
        α = line_search(policy_cache, data, θ, Δ)
        θ += α * Δ
        println("$i $(scn(policy_cache.c))  $(scn(α)) $(scn(norm(Δ,Inf)))")
    end

    return θ
end

function line_search(policy_cache, data, θ, Δ)
    c0 = policy_cache.c
    α = 1.0
    for i = 1:10
        evaluation!(policy_cache, data, θ + α * Δ)
        c_cand = policy_cache.c
        (c_cand < c0) && break
        α *= 0.5
    end
    return α
end

newton_solve(policy_cache, data)



pi
∂pi∂θ
∂pi∂xv
∂pi∂xvθ

c
g
H



################################################################################
# Simulation
################################################################################
# ## Initial conditions
q1 = nominal_configuration(RD.centroidal_quadruped)
v1 = zeros(RD.centroidal_quadruped.nq)

# ## Time
h = 0.02
timestep = h
T = 100

# ## Simulator
s = Simulator(RD.centroidal_quadruped, T, h=h)
s.ip.opts.r_tol = 1e-7
s.ip.opts.κ_tol = 1e-5
s.ip.opts.undercut = Inf
# ## Simulate
RD.simulate!(s, q1, v1)
# ## Visualize
RD.visualize!(vis, s)
set_light!(vis)
set_floor!(vis)







################################################################################
# Dynamics Model
################################################################################
dynamics_model = Simulator(RD.centroidal_quadruped, 1, h=h)
dynamics_model.ip.opts.r_tol = 1e-5
dynamics_model.ip.opts.κ_tol = 1e-2
dynamics_model.ip.opts.undercut = 5.0

nq = dynamics_model.model.nq
nx = 2nq
nu = dynamics_model.model.nu
nw = dynamics_model.model.nw
nu_infeasible = 6

################################################################################
# Gait design
################################################################################
file = JLD2.jldopen(joinpath(@__DIR__, "centroidal_quadruped_sol.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]

################################################################################
# iLQR
################################################################################
# ## initialization
# ## horizon
T = length(x_ref)
parameters = [[x_ref[i]; u_ref[min(i,T-1)]] for i=1:T]
x1 = deepcopy(x_ref[1])
xT = deepcopy(x_ref[T])

RD.set_robot!(vis, dynamics_model.model, x1)
RD.set_robot!(vis, dynamics_model.model, xT)

# ## (2-layer) multi-layer perceptron policy
l_input = nx
# l_input = nq + 3
l1 = 6
l2 = nu - 6
nθ = l1 * l_input + l2 * l1

function feature(x)
    xbody = x[1:3]
    xorient = x[4:6]
    xfoot1 = x[7:9]
    xfoot2 = x[10:12]
    xfoot3 = x[13:15]
    xfoot4 = x[16:18]

    xf = [1.0; xbody[2:3]; xorient; xfoot1 - xbody; xfoot2 - xbody; xfoot3 - xbody; xfoot4 - xbody; x[nq .+ (1:nq)]]
    # xf = [1.0; xbody[2:3]; xorient; xfoot1 - xbody; xfoot2 - xbody;
        # xfoot3 - xbody; xfoot4 - xbody; x[nq .+ (1:3)]]
    return xf
end

function policy(θ, x, w)
    shift = 0
    # input
    input = (feature(x) - feature(w[1:nx])) # policy independent of the x position

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


# ## model
h = timestep

function f1(y, x, u, w)
    @views u_ctrl = u[1:nu]
    @views x_di = x[1:nx]
    @views θ = u[nu .+ (1:nθ)]
    RD.dynamics(dynamics_model, view(y, 1:nx), x_di, u_ctrl, w)
    @views y[nx .+ (1:nθ)] .= θ
    return nothing
end

function f1x(dx, x, u, w)
    @views u_ctrl = u[1:nu]
    @views x_di = x[1:nx]
    @views θ = u[nu .+ (1:nθ)]
    dx .= 0.0
    RD.dynamics_jacobian_state(dynamics_model, view(dx, 1:nx, 1:nx), x_di, u_ctrl, w)
    return nothing
end

function f1u(du, x, u, w)
    @views u_ctrl = u[1:nu]
    @views x_di = x[1:nx]
    @views θ = u[nu .+ (1:nθ)]
    du .= 0.0
    RD.dynamics_jacobian_input(dynamics_model, view(du, 1:nx, 1:nu), x_di, u_ctrl, w)
    @views du[nx .+ (1:nθ), nu .+ (1:nθ)] .= I(nθ)
    return nothing
end

function ft(y, x, u, w)
    @views u_ctrl = u[1:nu]
    @views x_di = x[1:nx]
    @views θ = x[nx .+ (1:nθ)]
    RD.dynamics(dynamics_model, view(y, 1:nx), x_di, u_ctrl, w)
    @views y[nx .+ (1:nθ)] .= θ
    return nothing
end

function ftx(dx, x, u, w)
    @views u_ctrl = u[1:nu]
    @views x_di = x[1:nx]
    @views θ = x[nx .+ (1:nθ)]
    dx .= 0.0
    RD.dynamics_jacobian_state(dynamics_model, view(dx, 1:nx, 1:nx), x_di, u_ctrl, w)
    @views dx[nx .+ (1:nθ), nx .+ (1:nθ)] .= I(nθ)
    return nothing
end

function ftu(du, x, u, w)
    @views u_ctrl = u[1:nu]
    @views x_di = x[1:nx]
    @views θ = x[nx .+ (1:nθ)]
    RD.dynamics_jacobian_input(dynamics_model, view(du, 1:nx, 1:nu), x_di, u_ctrl, w)
    return nothing
end


# user-provided dynamics and gradients
dyn1 = iLQR.Dynamics(f1, f1x, f1u, nx + nθ, nx, nu + nθ)
dynt = iLQR.Dynamics(ft, ftx, ftu, nx + nθ, nx + nθ, nu)

dyn = [dyn1, [dynt for t = 2:T-1]...]

# ## objective
function o1(x, u, w)
    J = 0.0
    # qbody = [1e-0, 1e-0, 1e+1]
    # qfoot = [1e-0, 1e-0, 1e+2]
    # q = 1e-0 * [1e-0*qbody; 1e+0*ones(3); [qfoot; qfoot; qfoot; qfoot]]
    # v = 1e-0 * [1e-3*ones(3); 1e-2*ones(3); 1e-3*ones(12)]
    # r = 1e-2 * [ones(6); [1e-1,1,1e-2]; [1e-1,1,1e-2]; [1e-1,1,1e-2]; [1e-1,1,1e-2]]
    q = 1e+1 * ones(nq)
    v = 1e+1 * ones(nq)
    r = 1e+1 * ones(nu)
    ex = x[1:nx] - w[1:nx]
    eu = u[1:nu] - w[nx .+ (1:nu)]
    J += 0.5 * transpose(ex) * Diagonal([q; v]) * ex
    J += 0.5 * transpose(eu) * Diagonal(r) * eu
    J += 1e-1 * dot(u[nu .+ (1:nθ)], u[nu .+ (1:nθ)])
    return J
end

function ot(x, u, w)
    J = 0.0
    # qbody = [1e-0, 1e-0, 1e+1]
    # qfoot = [1e-0, 1e-0, 1e+2]
    # q = 1e-0 * [1e-0*qbody; 1e+0*ones(3); [qfoot; qfoot; qfoot; qfoot]]
    # v = 1e-0 * [1e-3*ones(3); 1e-2*ones(3); 1e-3*ones(12)]
    # r = 1e-2 * [ones(6); [1e-1,1,1e-2]; [1e-1,1,1e-2]; [1e-1,1,1e-2]; [1e-1,1,1e-2]]
    q = 1e+1 * ones(nq)
    v = 1e+1 * ones(nq)
    r = 1e+1 * ones(nu)
    ex = x[1:nx] - w[1:nx]
    eu = u[1:nu] - w[nx .+ (1:nu)]
    J += 0.5 * transpose(ex) * Diagonal([q; v]) * ex
    J += 0.5 * transpose(eu) * Diagonal(r) * eu
    J += 1e-1 * dot(x[nx .+ (1:nθ)], x[nx .+ (1:nθ)])
    return J
end

function oT(x, u, w)
    J = 0.0
    return J
end

c1 = iLQR.Cost(o1, nx, nu + nθ, num_parameter=nx+nu)
ct = iLQR.Cost(ot, nx + nθ, nu, num_parameter=nx+nu)
cT = iLQR.Cost(oT, nx + nθ, 0, num_parameter=nx+nu)
obj = [c1, [ct for t = 2:(T - 1)]..., cT]


# ## constraints
ul = -1.0 * [1e-1*ones(nu_infeasible); 1e3ones(nu-nu_infeasible)]
uu = +1.0 * [1e-1*ones(nu_infeasible); 1e3ones(nu-nu_infeasible)]

function con1(x, u, w)
    θ = u[nu .+ (1:nθ)]
    [
        1e+0 * (ul - u[1:nu]);
        1e+0 * (u[1:nu] - uu);
        1.0e-2 * (u[nu_infeasible+1:nu] - policy(θ, x[1:nx], w));
    ]
end

function cont(x, u, w)
    θ = x[nx .+ (1:nθ)]
    [
        1e+0 * (ul - u[1:nu]);
        1e+0 * (u[1:nu] - uu);
        1.0e-2 * (u[nu_infeasible+1:nu] - policy(θ, x[1:nx], w))
    ]
end

function goal(x, u, w)
    [
        # x[[1,nq+1]] - xT[[1,nq+1]];
        # x[nq+1:nq+1] - xT[nq+1:nq+1];
        1e+0 * (x[1:nq+1] - xT[1:nq+1]);
    ]
end

con_policy1 = iLQR.Constraint(con1, nx, nu + nθ, num_parameter=nx+nu, indices_inequality=collect(1:2nu))
con_policyt = iLQR.Constraint(cont, nx + nθ, nu, num_parameter=nx+nu, indices_inequality=collect(1:2nu))
con_policyT = iLQR.Constraint(goal, nx + nθ, 0)

cons = [con_policy1, [con_policyt for t = 2:T-1]..., con_policyT]
# ## problem
opts = iLQR.Options(line_search=:armijo,
    max_iterations=75,
    max_dual_updates=30,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-3,
    initial_constraint_penalty=1e-3,
    scaling_penalty=3.0,
    max_penalty=1e4,
    verbose=true)

p = iLQR.Solver(dyn, obj, cons, options=opts, parameters=parameters)

# ## initialize
θ0 = 1.0 * randn(nθ)
u_guess = [t == 1 ? [u_hover; θ0] : u_hover for t = 1:T-1]
x_guess = iLQR.rollout(dyn, x1, u_guess, parameters)

s = Simulator(RD.centroidal_quadruped, T-1, h=h)
for i = 1:T
    q = x_guess[i][1:nq]
    v = x_guess[i][nq .+ (1:nq)]
    RD.set_state!(s, q, v, i)
end
visualize!(vis, s)
# vis = Visualizer()
# open(vis)

iLQR.initialize_controls!(p, u_guess)
iLQR.initialize_states!(p, x_guess)
dynamics_model.ip.opts.r_tol = 1e-6
dynamics_model.ip.opts.κ_tol = 1e-2
local_continuation_callback!(solver::Solver) = continuation_callback!(solver, dynamics_model)

# ## solve
@time iLQR.constrained_ilqr_solve!(p, augmented_lagrangian_callback! = local_continuation_callback!)


# ## solution
x_sol, u_sol = iLQR.get_trajectory(p)
θ_sol = u_sol[1][nu .+ (1:nθ)]

# ## state
plot(hcat([x[1:nx] for x in x_sol]...)', label="", color=:orange, width=2.0)

# ## control
plot(hcat([u[1:nu] for u in u_sol]..., u_sol[end])', linetype = :steppost)

# ## plot xy
plot([x[1] for x in x_sol], [x[2] for x in x_sol], label="", color=:black, width=2.0)

# ## visualization
s = Simulator(RD.centroidal_quadruped, T-1, h=h)
for i = 1:T
    q = x_sol[i][1:nq]
    v = x_sol[i][nq .+ (1:nq)]
    RD.set_state!(s, q, v, i)
end
visualize!(vis, s)

# ## simulate policy
x_hist = [x1]
u_hist = [u_hover]

for t = 1:10T
    push!(u_hist, [zeros(nu_infeasible); policy(θ_sol, x_hist[end], gait[(t-1)%T+1])])
    y = zeros(nx)
    RD.dynamics(dynamics_model, y, x_hist[end], u_hist[end], zeros(nw))
    push!(x_hist, y)
end

s = Simulator(RD.centroidal_quadruped, 10T-1, h=h)
for i = 1:10T
    q = x_hist[i][1:nq]
    v = x_hist[i][nq .+ (1:nq)]
    RD.set_state!(s, q, v, i)
end
visualize!(vis, s)
set_light!(vis)
set_floor!(vis)

# Dojo.convert_frames_to_video_and_gif("RD.centroidal_quadruped_single_regularized_open_loop")
# Dojo.convert_frames_to_video_and_gif("RD.centroidal_quadruped_single_regularized_policy")
