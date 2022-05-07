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
# Simple policies
################################################################################

function pd_policy(x, θ, x_ref, u_ref)
    q = x[1:nq]
    Δx = x_ref - x
    J = input_jacobian(dynamics_model.model, q)
    Δxq = Δx[1:nq]
    Δxv = Δx[nq .+ (1:nq)]
    Kp = 1e-0
    Kv = 1e-1
    u = u_ref + Kp * (J \ Δxq) + Kv * (J \ Δxv)
    return u
end

function K_policy(x, θ, x_ref, u_ref, K_ref)
    q = x[1:nq]
    Δx = x - x_ref
    u = u_ref + K_ref * Δx
    u[1:6] .= 0.0
    return u
end

file = JLD2.jldopen(joinpath(@__DIR__, "dataset/centroidal_quadruped_ref.jld2"))
x_ref = file["x_sol"]
u_ref = file["u_sol"]
K_ref = file["K_sol"]
close(file)

# stride
Δx = x_ref[end-1] - x_ref[1]
x_stride = mean(Δx[[1,7,10,13,16]])

γ = 0.3
x_hist = [x_ref[1] + [0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; 0;0;γ; zeros(nq)]]
u_hist = []

M = 20
t = 0
for i = 1:M
    for j = 1:T-1
        t += 1
        t_loop = (t-1)%(T-1)+1
        x_offset = deepcopy(x_ref[t_loop])
        x_offset[[1,7,10,13,16]] .+= (i-1) * x_stride
        x_offset[1 .+ [1,7,10,13,16]] .+= (i-1) * 0.01
        push!(u_hist, K_policy(x_hist[end], θ0, x_offset, u_ref[t_loop], K_ref[t_loop]))
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
