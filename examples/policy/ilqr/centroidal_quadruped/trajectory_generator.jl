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

number_sample = 40
trajectory_sampler(ilqr_solver, dynamics_model, gait, number_sample;
        u_guess=[[4,0,16, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0.0] for i=1:length(gait)-1],
        configuration_amplitude=0.10,
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
    # plt = plot(hcat(u_sols[i]...)'[:,1:nu], color=:red, legend=false)
    # display(plt)
    # sleep(0.1)
end
plot!(plt, hcat(u_ref...)'[:,1:nu], linewidth=3.0, color=:black, legend=false)
display(plt)

plt = plot()
for i = 1:number_sample
    plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_sols[i]]...)'[:,1:nu], color=:red, legend=false)
end
plot!(plt, hcat([reshape(K, (nu*nx)) for K in K_ref]...)'[:,1:nu], linewidth=3.0, color=:black, legend=false)
display(plt)


[x[end] for x in x_sols]

anim = MeshCat.Animation(convert(Int, floor(1.0 / 0.02)))
for (i,x_sol) in enumerate(x_sols)
    s = Simulator(RD.centroidal_quadruped, T-1, h=h)
    anim = visualize!(vis, s.model, [x[1:nq] for x in x_sol], name=Symbol(:robot, i), anim=anim)
    set_light!(vis)
    set_floor!(vis)
end
