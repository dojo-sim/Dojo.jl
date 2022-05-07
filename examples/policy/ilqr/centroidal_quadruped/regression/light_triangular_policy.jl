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
include("../../robodojo/codegen.jl")
RoboDojo.RESIDUAL_EXPR


include("../utils.jl")
include("../gait_design.jl")



################################################################################
# Time-varying Linear Least Squares
################################################################################

L = 1

function linear_policy(x, θ, x_ref, u_ref, i)
    ii = min(i, L)
    θmat = reshape(θ, (nu, ii*nx))
    Δx = x - x_ref
    u = u_ref + θmat * Δx
    return u
end

function linear_policy_jacobian_state(θ, i)
    ii = min(i, L)
    θmat = reshape(θ, (nu, ii*nx))
    return θmat
end

function linear_policy_jacobian_last_state(θ)
    θmat = reshape(θ[end-nu*nx+1:end], (nu, nx))
    return θmat
end

∂policy∂xθ0 = sparse(FiniteDiff.finite_difference_jacobian(θ -> reshape(linear_policy_jacobian_last_state(θ, 3), nu*nx), rand(nu*3nx)))
∂policy∂xθ = [spzeros(nu*nx, nu*(3-1)*nx) I(nu*1nx)]
function linear_policy_jacobian_last_state_parameters(i)
    ii = min(i, L)
    return [spzeros(nu*nx, nu*(ii-1)*nx) I(nu*1nx)]
end

function linear_policy_jacobian_parameters(x, x_ref, i)
    ii = min(i, L)
    Δx = x - x_ref
    ∂policy∂θ = hcat([Δx[j] * I(nu) for j=1:ii*nx]...)
    return ∂policy∂θ
end

# x0 = ones(2nx)
# J0 = linear_policy_jacobian_parameters(x0, vcat(x_ref[1:2]...), 2)
# J1 = FiniteDiff.finite_difference_jacobian(θ -> linear_policy(x0, θ, vcat(x_ref[1:2]...), u_ref[1], 2), rand(nu*2nx))
# norm(J0 - J1)

function evaluation(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols;
        Qu=Qu0*I(nu)/nu, Qθ=1e-10, QK=QK0*I(nu*nx)/nu/nx)

    sample_number = length(x_sols)
    ii = min(i, L)
    c = 0.0

    for j = 1:sample_number
        x_sol = vcat(x_sols[j][i-ii+1:i]...)
        u_sol = u_sols[j][i]
        K_sol = K_sols[j][i]
        # x_sol = x_ref[i]
        # u_sol = u_ref[i]
        # K_sol = K_ref[i]
        xi = vcat(x_ref[i-ii+1:i]...)
        ui = u_ref[i]
        Δu = linear_policy(x_sol, θ, xi, ui, i) - u_sol
        # ΔK = reshape(linear_policy_jacobian_last_state(θ) - K_ref[i], nu*nx)
        ΔK = reshape(linear_policy_jacobian_last_state(θ) - K_sol, nu*nx)
        c += 0.5 * Δu' * Qu * Δu
        c += 0.5 * ΔK' * QK * ΔK
    end
    c /= sample_number
    c += 0.5 * θ' * Qθ * θ
    return c
end

function gradient(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols;
        Qu=Qu0*I(nu)/nu, Qθ=1e-10, QK=QK0*I(nu*nx)/nu/nx)

    sample_number = length(x_sols)
    ii = min(i, L)
    grad = zeros(nu*ii*nx)

    for j = 1:sample_number
        x_sol = vcat(x_sols[j][i-ii+1:i]...)
        u_sol = u_sols[j][i]
        K_sol = K_sols[j][i]
        # x_sol = x_ref[i]
        # u_sol = u_ref[i]
        # K_sol = K_ref[i]
        xi = vcat(x_ref[i-ii+1:i]...)
        ui = u_ref[i]
        Δu = linear_policy(x_sol, θ, xi, ui, i) - u_sol
        Δu∂θ = linear_policy_jacobian_parameters(x_sol, xi, i)
        # ΔK = reshape(linear_policy_jacobian_last_state(θ) - K_ref[i], nu*nx)
        ΔK = reshape(linear_policy_jacobian_last_state(θ) - K_sol, nu*nx)
        ΔK∂θ = linear_policy_jacobian_last_state_parameters(i)
        grad += Δu∂θ' * Qu * Δu
        grad += ΔK∂θ' * QK * ΔK
    end
    grad /= sample_number
    grad += Qθ * θ
    return grad
end

function hessian(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols;
        Qu=Qu0*I(nu)/nu, Qθ=1e-10, QK=QK0*I(nu*nx)/nu/nx)

    sample_number = length(x_sols)
    ii = min(i, L)
    hess = spzeros(nu*ii*nx, nu*ii*nx)

    for j = 1:sample_number
        x_sol = vcat(x_sols[j][i-ii+1:i]...)
        u_sol = u_sols[j][i]
        K_sol = K_sols[j][i]
        # x_sol = x_ref[i]
        # u_sol = u_ref[i]
        # K_sol = K_ref[i]
        xi = vcat(x_ref[i-ii+1:i]...)
        ui = u_ref[i]
        Δu∂θ = linear_policy_jacobian_parameters(x_sol, xi, i)
        ΔK∂θ = linear_policy_jacobian_last_state_parameters(i)
        hess += Δu∂θ' * Qu * Δu∂θ
        hess += ΔK∂θ' * QK * ΔK∂θ
    end
    hess /= sample_number
    hess += Qθ * I
    return hess
end


number_sample = 40
x_sols = []
u_sols = []
K_sols = []
for i = 1:number_sample
    file = JLD2.jldopen(joinpath(@__DIR__, "../dataset/centroidal_quadruped_$i.jld2"))
    push!(x_sols, file["x_sol"])
    push!(u_sols, file["u_sol"])
    push!(K_sols, file["K_sol"])
    close(file)
end
file = JLD2.jldopen(joinpath(@__DIR__, "../dataset/centroidal_quadruped_ref.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]
K_ref = file["K_sol"]
close(file)

i0 = 11
ii0 = min(i0, L)
θ0 = zeros(nu*ii0*nx)
evaluation(θ0, i0, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
G0 = gradient(θ0, i0, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
H0 = hessian(θ0, i0, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
θ0 = θ0 - H0 \ G0
evaluation(θ0, i0, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)


Qu0 = 1e-10
QK0 = 1e-1
K_fit = [zeros(nu,min(i,L)*nx) for i=1:T-1]
for i = 1:T-1
    ii = min(i, L)
    θ = zeros(nu*ii*nx)
    evaluation(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
    G = gradient(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
    H = hessian(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
    θ = θ - H \ G
    E = evaluation(θ, i, x_ref, u_ref, K_ref, x_sols, u_sols, K_sols)
    @show i, E
    K_fit[i] .= reshape(θ, (nu,ii*nx))
end



file = JLD2.jldopen(joinpath(@__DIR__, "../dataset/centroidal_quadruped_ref.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]
K_ref = file["K_sol"]
close(file)

# stride
Δx = x_ref[end-1] - x_ref[1]
x_stride = mean(Δx[[1,7,10,13,16]])

γ = 0.01
x_hist = [x_ref[1] + [0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; zeros(nq)]]
u_hist = []

M = 10
t = 0
for i = 1:M
    for j = 1:T-1
        jj = min(j, L)
        t += 1
        t_loop = (t-1)%(T-1)+1
        x_offset = vcat(deepcopy(x_ref[t_loop-jj+1:t_loop])...)
        x_offset[[1,7,10,13,16]] .+= (i-1) * x_stride
        x_offset[1 .+ [1,7,10,13,16]] .+= (i-1) * 0.00
        # push!(u_hist, linear_policy(x_hist[end], K_sols[2][j], x_offset, u_ref[t_loop]))
        # push!(u_hist, [0;0;-0; zeros(3); linear_policy(x_hist[end], K_fit[j], x_offset, u_ref[t_loop])[7:end]])
        push!(u_hist, [0;0;-0; zeros(3); linear_policy(x_hist[end], K_ref[j], x_offset, u_ref[t_loop])[7:end]])
        # push!(u_hist, linear_policy(vcat(x_hist[end-jj+1:end]...), K_fit[j], x_offset, u_ref[t_loop], j))
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
    plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_sols[i]]...)'[:,1:nu],
        color=:red, legend=false)
end
plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_ref]...)'[:,1:nu],
    linewidth=3.0, color=:black, legend=false)
plot!(plt, hcat([reshape(K_fit[i][end-nu*nx+1:end], (nu*nx)) for i = 1:T-1]...)'[:,1:nu],
    linewidth=3.0, color=:lightblue, legend=false)
display(plt)
