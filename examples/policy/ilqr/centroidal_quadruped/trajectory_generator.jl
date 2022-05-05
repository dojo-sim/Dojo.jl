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


include("utils.jl")
include("gait_design.jl")

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
reset!(dynamics_model)
nq = dynamics_model.model.nq
nx = 2nq
nu = dynamics_model.model.nu
nw = dynamics_model.model.nw
nu_infeasible = 6

################################################################################
# Gait design
################################################################################
T = Int(floor(0.65 / h)) + 1
Tm = Int((T + 1) / 2)
s = Simulator(RD.centroidal_quadruped, T, h=h)

gait = trotting_gait(centroidal_quadruped, Tm, timestep=timestep, velocity=0.15)
for x in gait
    RD.set_robot!(vis, centroidal_quadruped, x[1:nq])
    sleep(h)
end


################################################################################
# iLQR sampler
################################################################################

ilqr_solver = solver_generator(dynamics_model, gait;
    u_hover=[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0])
x_sol, u_sol, K_sol, max_violation = trajectory_generator(ilqr_solver, dynamics_model, gait,
    u_guess=[[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0] for i=1:length(gait)-1],
    initial_disturbance=[+0.0;0;0.0; 0;0;0; +0.0;0;0.0; +0.0;0;0.0; +0.0;0;0.0; +0.0;0;0.0; zeros(nq)]
    )

JLD2.jldsave(joinpath(@__DIR__, "dataset/centroidal_quadruped_ref.jld2"),
    x_sol=x_sol, u_sol=u_sol, K_sol=K_sol)

sample_number = 40
trajectory_sampler(ilqr_solver, dynamics_model, gait, sample_number;
        u_guess=[[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0] for i=1:length(gait)-1],
        configuration_amplitude=0.15,
        velocity_amplitude=0.00,
        )



x_sols = []
u_sols = []
K_sols = []
for i = 1:number_sample
    file = JLD2.jldopen(joinpath(@__DIR__, "dataset/centroidal_quadruped_$i.jld2"))
    push!(x_sols, file["x_sol"])
    push!(u_sols, file["u_sol"])
    push!(K_sols, file["K_sol"])
    close(file)
end
file = JLD2.jldopen(joinpath(@__DIR__, "dataset/centroidal_quadruped_ref.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]
K_ref = file["K_sol"]
close(file)

plt = plot()
for i = 1:number_sample
    plot!(plt, hcat(x_sols[i]...)'[:,1:nq], color=:red, legend=false)
end
plot!(plt, hcat(x_ref...)'[:,1:nq], linewidth=3.0, color=:black, legend=false)
display(plt)

plt = plot()
for i = 1:number_sample
    plot!(plt, hcat(u_sols[i]...)'[:,1:nu], color=:red, legend=false)
end
plot!(plt, hcat(u_ref...)'[:,1:nu], linewidth=3.0, color=:black, legend=false)
display(plt)
display(plt)

################################################################################
# Policy
################################################################################

function feature(x)
    nq = 18
    xbody = x[1:3]
    xorient = x[4:6]
    xfoot1 = x[7:9]
    xfoot2 = x[10:12]
    xfoot3 = x[13:15]
    xfoot4 = x[16:18]

    # xf = [1.0; xbody[1:3]; xorient;
    #     xfoot1 - xbody;
    #     xfoot2 - xbody;
    #     xfoot3 - xbody;
    #     xfoot4 - xbody;
    #     # x[nq .+ (1:nq)]
    #     ]
    xf = x[1:nq]
    return xf
end

nf = length(feature(zeros(nx)))
∂feature∂x = FiniteDiff.finite_difference_jacobian(x -> feature(x), rand(nx))
function feature_jacobian_state()
    return ∂feature∂x
end

function policy(x, θ, x_ref, u_ref)
    θmat = reshape(θ, (nu,nf))
    u = u_ref + θmat * (feature(x) - feature(x_ref))
    return u
end

function policy_jacobian_state(θ)
    θmat = reshape(θ, (nu,nf))
    return θmat * feature_jacobian_state()
end

∂policy∂xθ = FiniteDiff.finite_difference_jacobian(θ -> reshape(policy_jacobian_state(θ), nu*nx), rand(nu*nf))
function policy_jacobian_state_parameters()
    return ∂policy∂xθ
end

function policy_jacobian_parameters(x, x_ref)
    f = feature(x) - feature(x_ref)
    ∂policy∂θ = hcat([f[i] * I(nu) for i=1:nf]...)
    return ∂policy∂θ
end

u_ref0 = ones(nu)
x_ref0 = ones(nx)
x0 = ones(nx)
θ0 = rand(nu*nf)
policy(x0, θ0, x_ref0, u_ref0)
policy_jacobian_parameters(x0, x_ref0)
policy_jacobian_state(θ0)

function evaluation(θ, x_ref, u_ref, x_sols, u_sols, K_sols;
        Qu=I(nu), Qθ=1e-8*I(nu*nf), QK=1e-7*I(nu*nx)/nu/nx, sample_number=20, batch=5)
    T = length(x_ref)
    counter = 0
    c = 0.0

    for i = 1:sample_number#in rand(1:sample_number, batch)
        x_sol = x_sols[i]
        u_sol = u_sols[i]
        K_sol = K_sols[i]
        for i = 1:T-1
            Δu = policy(x_sol[i], θ, x_ref[i], u_ref[i]) - u_sol[i]
            ΔK = reshape(policy_jacobian_state(θ) - K_sol[i], nu*nx)
            c += 0.5 * Δu' * Qu * Δu
            c += 0.5 * ΔK' * QK * ΔK
            counter += 1
        end
    end
    c /= counter
    c += 0.5 * θ' * Qθ * θ
    return c
end

function gradient(θ, x_ref, u_ref, x_sols, u_sols, K_sols;
        Qu=I(nu), Qθ=1e-8*I(nu*nf), QK=1e-7*I(nu*nx)/nu/nx, sample_number=20, batch=5)
    T = length(x_ref)
    counter = 0
    grad = zeros(nu*nf)

    for i = 1:sample_number#in rand(1:sample_number, batch)
        x_sol = x_sols[i]
        u_sol = u_sols[i]
        K_sol = K_sols[i]
        for i = 1:T-1
            Δu = policy(x_sol[i], θ, x_ref[i], u_ref[i]) - u_sol[i]
            Δu∂θ = policy_jacobian_parameters(x_sol[i], x_ref[i])
            ΔK = reshape(policy_jacobian_state(θ) - K_sol[i], nu*nx)
            ΔK∂θ = policy_jacobian_state_parameters()
            grad += Δu∂θ' * Qu * Δu
            grad += ΔK∂θ' * QK * ΔK
            counter += 1
        end
    end
    grad /= counter
    grad += Qθ * θ
    return grad
end

function hessian(θ, x_ref, u_ref, x_sols, u_sols, K_sols;
        Qu=I(nu), Qθ=1e-8*I(nu*nf), QK=1e-7*I(nu*nx)/nu/nx, sample_number=20, batch=5)
    T = length(x_ref)
    counter = 0
    hess = zeros(nu*nf, nu*nf)

    for i = 1:sample_number#in rand(1:sample_number, batch)
        x_sol = x_sols[i]
        u_sol = u_sols[i]
        K_sol = K_sols[i]
        for i = 1:T-1
            Δu∂θ = policy_jacobian_parameters(x_sol[i], x_ref[i])
            ΔK∂θ = policy_jacobian_state_parameters()
            hess += Δu∂θ' * Qu * Δu∂θ
            hess += ΔK∂θ' * QK * ΔK∂θ
            counter += 1
        end
    end
    hess /= counter
    hess += Qθ
    return hess
end


function gradient_descent(θ, x_ref, u_ref, x_sols, u_sols, K_sols;
        number_iterations=10, α=1.0, batch=5)
    eval_prev = -Inf
    θ_prev = deepcopy(θ)

    hess = hessian(θ, x_ref, u_ref, x_sols, u_sols, K_sols, batch=batch) # constant

    for i = 1:number_iterations
        grad = gradient(θ, x_ref, u_ref, x_sols, u_sols, K_sols, batch=batch)
        θ = θ - α * (hess+ 1e-6 * I) \ grad
        if (i-1) % 1 == 0
            eval = evaluation(θ, x_ref, u_ref, x_sols, u_sols, K_sols, batch=batch)
            @show i, eval, α
            if eval_prev > eval
                α = clamp(α * 2.0, 1e-4, 1e0)
            else
                α = clamp(α / 1.5, 1e-4, 1e0)
                θ = θ_prev
            end
            eval_prev = eval
            θ_prev .= θ
        end
    end
    return θ
end

number_sample = 20
x_sols = []
u_sols = []
K_sols = []
for i = 1:number_sample
    file = JLD2.jldopen(joinpath(@__DIR__, "dataset/centroidal_quadruped_$i.jld2"))
    push!(x_sols, file["x_sol"])
    push!(u_sols, file["u_sol"])
    push!(K_sols, file["K_sol"])
    close(file)
end
file = JLD2.jldopen(joinpath(@__DIR__, "dataset/centroidal_quadruped_ref.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]
K_ref = file["K_sol"]
close(file)

# θ0 = -3e-2 * rand(nu*nf)
θ0 = gradient_descent(θ0, x_ref, u_ref, x_sols, u_sols, K_sols,
    number_iterations=100, batch=20, α=1.0)
θ1 = deepcopy(θ0)


# θ0 = zeros(nu*nf)
# ## simulate policy

γ = 0.00
x_hist = [x_ref[1] + [0;0;γ; 0;0;0; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; zeros(nq)]]
u_hist = []

for t = 1:3T-1
    t_loop = (t-1)%(T-1)+1
    # push!(u_hist, [zeros(nu_infeasible); policy(θ_sol, x_hist[end], gait[(t-1)%T+1])])
    push!(u_hist, policy(x_hist[end], θ0, x_ref[t_loop], u_ref[t_loop]))
    # @show policy(x_hist[end], θ0, x_ref[t_loop], u_ref[t_loop]) - u_ref[t_loop]
    # push!(u_hist, u_ref[t_loop])
    y = zeros(nx)
    RD.dynamics(dynamics_model, y, x_hist[end], u_hist[end], zeros(nw))
    push!(x_hist, y)
end

s = Simulator(RD.centroidal_quadruped, 3T-1, h=h)
for i = 1:3T
    q = x_hist[i][1:nq]
    v = x_hist[i][nq .+ (1:nq)]
    RD.set_state!(s, q, v, i)
end
visualize!(vis, s)
set_light!(vis)
set_floor!(vis)

plot(Gray.(abs.(reshape(θ0 ./ maximum(θ0), (nu,nf)))))
spy(reshape(θ0, (nu,nf)), markersize=5.0)




################################################################################
# LBFGS
################################################################################

using Optim

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

function g!(G, x)
    G[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
    G[2] = 200.0 * (x[2] - x[1]^2)
end

optimize(f, g!, x0, LBFGS())



function evaluation(θ, x_ref, u_ref, x_sols, u_sols, K_sols;
        Qu=I(nu), Qθ=1e-8*I(nu*nf), QK=1e-2*I(nu*nx)/nu/nx, sample_number=20, batch=5)
    T = length(x_ref)
    counter = 0
    c = 0.0

    for i = 1:sample_number#in rand(1:sample_number, batch)
        x_sol = x_sols[i]
        u_sol = u_sols[i]
        K_sol = K_sols[i]
        for i = 1:T-1
            Δu = policy(x_sol[i], θ, x_ref[i], u_ref[i]) - u_sol[i]
            ΔK = reshape(policy_jacobian_state(θ) - K_sol[i], nu*nx)
            c += 0.5 * Δu' * Qu * Δu
            c += 0.5 * ΔK' * QK * ΔK
            counter += 1
        end
    end
    c /= counter
    c += 0.5 * θ' * Qθ * θ
    return c
end

function gradient(θ, x_ref, u_ref, x_sols, u_sols, K_sols;
        Qu=I(nu), Qθ=1e-8*I(nu*nf), QK=1e-2*I(nu*nx)/nu/nx, sample_number=20, batch=5)
    T = length(x_ref)
    counter = 0
    grad = zeros(nu*nf)

    for i = 1:sample_number#in rand(1:sample_number, batch)
        x_sol = x_sols[i]
        u_sol = u_sols[i]
        K_sol = K_sols[i]
        for i = 1:T-1
            Δu = policy(x_sol[i], θ, x_ref[i], u_ref[i]) - u_sol[i]
            Δu∂θ = policy_jacobian_parameters(x_sol[i], x_ref[i])
            ΔK = reshape(policy_jacobian_state(θ) - K_sol[i], nu*nx)
            ΔK∂θ = policy_jacobian_state_parameters()
            grad += Δu∂θ' * Qu * Δu
            grad += ΔK∂θ' * QK * ΔK
            counter += 1
        end
    end
    grad /= counter
    grad += Qθ * θ
    return grad
end
