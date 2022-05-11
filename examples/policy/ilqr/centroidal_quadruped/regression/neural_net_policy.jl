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

include("../../methods.jl")

vis = Visualizer()
open(vis)

include("../gait_design.jl")

include("../../../robodojo/centroidal_quadruped/model.jl")
include("../../../robodojo/centroidal_quadruped/visuals.jl")
include("../../../robodojo/centroidal_quadruped/simulator.jl")
include("../../../robodojo/dynamics.jl")

RoboDojo.RESIDUAL_EXPR
force_codegen = true
# force_codegen = false
robot = centroidal_quadruped
include("../../../robodojo/codegen.jl")
RoboDojo.RESIDUAL_EXPR


include("../utils.jl")
include("../gait_design.jl")


################################################################################
# Time-varying Linear Least Squares
################################################################################

using NablaNet


function policy(x, θ, x_ref, u_ref, K_ref, net)
    Δx = x - x_ref
    NablaNet.evaluation!(net, Δx, θ)
    u_net = NablaNet.get_output(net)
    u = u_ref + K_ref * Δx + u_net
    return u
end

function policy_jacobian_state(x, θ, x_ref, u_ref, K_ref, net)
    Δx = x - x_ref
    NablaNet.jacobian_input!(net, Δx, θ)
    ∂u_net∂x = net.jacobian_input
    Jx = K_ref + ∂u_net∂x
    return Jx
end

function policy_jacobian_parameters(x, θ, x_ref, u_ref, K_ref, net)
    Δx = x - x_ref
    NablaNet.jacobian_parameters!(net, Δx, θ)
    ∂u_net∂θ = net.jacobian_parameters
    Jθ = ∂u_net∂θ
    return Jθ
end

x0 = rand(nx)
x_ref0 = rand(nx)
u_ref0 = rand(nu)
K_ref0 = rand(nu,nx)
# net0 = NablaNet.Net(nx, nu, dim_layers=[30,20,18], activations=[x->log.(1 .+ exp.(x)), x->log.(1 .+ exp.(x)), x->log.(1 .+ exp.(x)), x->x])
# net0 = NablaNet.Net(nx, nu, dim_layers=[20], activations=[x->log.(1 .+ exp.(x)), x->x])
net0 = NablaNet.Net(nx, nu, dim_layers=[20], activations=[x->tanh.(x), x->x])
# net0 = NablaNet.Net(nx, nu)
nθ = NablaNet.parameter_dimension(net0)
θ0 = rand(nθ)
u0 = policy(x0, θ0, x_ref0, u_ref0, K_ref0, net0)
# @benchmark policy(x0, θ0, x_ref0, u_ref0, K_ref0, net0)

Jx1 = policy_jacobian_state(x0, θ0, x_ref0, u_ref0, K_ref0, net0)
# @benchmark policy_jacobian_state(x0, θ0, x_ref0, u_ref0, K_ref0, net0)
Jx0 = FiniteDiff.finite_difference_jacobian(x0 -> policy(x0, θ0, x_ref0, u_ref0, K_ref0, net0), x0)
norm(Jx1 - Jx0, Inf)

Jθ1 = policy_jacobian_parameters(x0, θ0, x_ref0, u_ref0, K_ref0, net0)
# @benchmark policy_jacobian_parameters(x0, θ0, x_ref0, u_ref0, K_ref0, net0)
Jθ0 = FiniteDiff.finite_difference_jacobian(θ0 -> policy(x0, θ0, x_ref0, u_ref0, K_ref0, net0), θ0)
norm(Jθ1 - Jθ0, Inf)


function evaluation(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols, net;
        Qu=1e-0, Qθ=1e-3, QK=1e-0)

    sample_number = length(x_sols)
    c = 0.0

    for j = 1:sample_number
        x_sol = x_sols[j][i]
        u_sol = u_sols[j][i]
        K_sol = K_sols[j][i]
        xi = x_ref[i]
        ui = u_ref[i]
        Ki = K_ref[i]
        Δu = policy(x_sol, θ, xi, ui, Ki, net) - u_sol
        # ΔK = reshape(linear_policy_jacobian_state(θ) - K_sol, nu*nx)
        c += 0.5 * Δu' * Qu * Δu
        # c += 0.5 * ΔK' * QK * ΔK
    end
    c /= sample_number
    c += 0.5 * θ' * Qθ * θ
    return c
end

function gradient(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols, net;
        Qu=1e-0, Qθ=1e-3, QK=1e-0)

    sample_number = length(x_sols)
    nθ = NablaNet.parameter_dimension(net)
    grad = zeros(nθ)

    for j = 1:sample_number
        x_sol = x_sols[j][i]
        u_sol = u_sols[j][i]
        K_sol = K_sols[j][i]
        xi = x_ref[i]
        ui = u_ref[i]
        Ki = K_ref[i]
        Δu = policy(x_sol, θ, xi, ui, Ki, net) - u_sol
        Δu∂θ = policy_jacobian_parameters(x_sol, θ, xi, ui, Ki, net)
        # ΔK = reshape(linear_policy_jacobian_state(θ) - K_sol, nu*nx)
        # ΔK∂θ = linear_policy_jacobian_state_parameters()
        grad += Δu∂θ' * Qu * Δu
        # grad += ΔK∂θ' * QK * ΔK
    end
    grad /= sample_number
    grad += Qθ * θ
    return grad
end

function hessian(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols, net;
        Qu=1e-0, Qθ=1e-3, QK=1e-0)

    sample_number = length(x_sols)
    nθ = NablaNet.parameter_dimension(net)
    hess = zeros(nθ, nθ)

    for j = 1:sample_number
        x_sol = x_sols[j][i]
        u_sol = u_sols[j][i]
        K_sol = K_sols[j][i]
        xi = x_ref[i]
        ui = u_ref[i]
        Ki = K_ref[i]
        # Δu∂θ = linear_policy_jacobian_parameters(x_sol, xi)
        Δu∂θ = policy_jacobian_parameters(x_sol, θ, xi, ui, Ki, net)
        # ΔK∂θ = linear_policy_jacobian_state_parameters()
        hess += Δu∂θ' * Qu * Δu∂θ
        # hess += ΔK∂θ' * QK * ΔK∂θ
    end
    hess /= sample_number
    hess += Qθ * I
    return hess
end


number_sample0 = 40
x_sols0 = []
u_sols0 = []
K_sols0 = []
for i = 1:number_sample0
    file = JLD2.jldopen(joinpath(@__DIR__, "../dataset/centroidal_quadruped_$i.jld2"))
    push!(x_sols0, file["x_sol"])
    push!(u_sols0, file["u_sol"])
    push!(K_sols0, file["K_sol"])
    close(file)
end
file = JLD2.jldopen(joinpath(@__DIR__, "../dataset/centroidal_quadruped_ref.jld2"))
x_ref0 = file["x_sol"]
u_ref0 = file["u_sol"]
K_ref0 = file["K_sol"]
close(file)

# for i = 1:number_sample
#     for j = 1:length(u_ref)
#         Δx = x_sols[i][j] - x_ref[j]
#         u_sols[i][j] = u_ref[j] + K_ref[j] * Δx
#     end
# end



# θ0 = reshape([zeros(nu, 1) deepcopy(K_ref[1])], nθ)
# θ0 += 1e-0 * (rand(nθ) .- 0.5)
θ0 = 1e-3 * (rand(nθ) .- 0.5)
i0 = 1
Qu0 = 1*1e-0
Qθ0 = 1*1e-8
evaluation(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0)
G0 = gradient(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0)
G1 = FiniteDiff.finite_difference_gradient(θ0 ->
    evaluation(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0), θ0)
norm(G0 - G1)


H0 = hessian(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0)
for j = 1:500
    E = evaluation(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0)
    Eu = evaluation(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=0*Qθ0)
    Eθ = evaluation(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=0*Qu0, Qθ=Qθ0)

    G0 = gradient(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0)
    (j % 25 == 0) && (H0 = hessian(θ0, i0, x_ref0, u_ref0, K_ref0, x_sols0, u_sols0, K_sols0, net0, Qu=Qu0, Qθ=Qθ0))
    θ0 = θ0 - 1e-2 * ((H0 + (1e-5+E)*1e-2 * I) \ G0)
    # θ0 = θ0 - 1e-2 * G0
    println("j ", j, "  E", scn(E), "  Eu", scn(Eu), "  Eθ", scn(Eθ), "  |G|", scn(norm(G0)), "  |θ|", scn(norm(θ0)))
end



u_sols0
Δu0 = policy(x_sols0[1][1], θ0, x_ref0[1], u_ref0[1], K_ref0[1], net0) - u_sols0[1][1]
norm(Δu0)
Δu0 = policy(x_sols0[2][1], θ0, x_ref0[1], u_ref0[1], K_ref0[1], net0) - u_sols0[2][1]
norm(Δu0)
Δu0 = policy(x_sols0[3][1], θ0, x_ref0[1], u_ref0[1], K_ref0[1], net0) - u_sols0[3][1]
norm(Δu0)


policy(x_sols0[1][1], θ0, x_ref0[1], u_ref0[1], K_ref0[1], net0)
policy(x_sols0[1][1] - ones(nx), θ0, x_ref0[1] - ones(nx), 1 .+ u_ref0[1], K_ref0[1], net0)
policy(x_sols0[2][1], θ0, x_ref0[1], u_ref0[1], K_ref0[1], net0)


a = 2
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1
a = 1












file = JLD2.jldopen(joinpath(@__DIR__, "../dataset/centroidal_quadruped_ref.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]
K_ref = file["K_sol"]
close(file)

# stride
Δx = x_ref[end-1] - x_ref[1]
x_stride = mean(Δx[[1,7,10,13,16]])

γ = 0.10
x_hist = [x_ref[1] + [0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; zeros(nq)]]
u_hist = []

M = 10
t = 0
for i = 1:M
    for j = 1:T-1
        t += 1
        t_loop = (t-1)%(T-1)+1
        x_offset = deepcopy(x_ref[t_loop])
        x_offset[[1,7,10,13,16]] .+= (i-1) * x_stride
        x_offset[1 .+ [1,7,10,13,16]] .+= (i-1) * 0.00
        # push!(u_hist, linear_policy(x_hist[end], K_sols[2][j], x_offset, u_ref[t_loop]))
        # push!(u_hist, [0;0;-0; zeros(3); linear_policy(x_hist[end], K_fit[j], x_offset, u_ref[t_loop])[7:end]])
        push!(u_hist, [0;0;-0; zeros(3); linear_policy(x_hist[end], K_ref[j], x_offset, u_ref[t_loop])[7:end]])
        # push!(u_hist, linear_policy(x_hist[end], K_fit[j], x_offset, u_ref[t_loop]))
        # push!(u_hist, linear_policy(x_hist[end], K_ref[j], x_offset, u_ref[t_loop]))
        y = zeros(nx)
        RD.dynamics(dynamics_model, y, x_hist[end], u_hist[end], zeros(nw))
        push!(x_hist, y)
    end
end

s = Simulator(RD.centroidal_quadruped, M*(T-1), h=h)
for i = 1:M*(T-1)+1
    q = x_hist[i][1:nq]
    v = x_hist[i][nq .+ (1:nq)]
    RD.set_state!(s, q, v, i)
end
visualize!(vis, s)
set_light!(vis)
set_floor!(vis)

plot(hcat(u_hist...)')

norm.(K_fit)
norm.(K_ref)
norm.(K_fit .- K_ref)
# norm.(2K_fit .- K_ref)

plt = plot()
for i = 1:number_sample
    plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_sols[i]]...)'[:,1:nu], color=:red, legend=false)
end
plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_ref]...)'[:,1:nu], linewidth=3.0, color=:black, legend=false)
plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_fit]...)'[:,1:nu], linewidth=3.0, color=:lightblue, legend=false)
display(plt)
