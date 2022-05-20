using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
Pkg.instantiate()

# ## visualizer
vis = Visualizer()
open(vis)

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra
using FiniteDiff
using DojoEnvironments
using JLD2

# ## scripts
include(joinpath(module_dir(), "examples/policy/methods/continuation.jl"))
include(joinpath(module_dir(), "examples/policy/methods/tvlqr.jl"))
include(joinpath(module_dir(), "DojoEnvironments/src",
    "quadruped/methods/template.jl"))


################################################################################
# ## system
################################################################################
gravity = -9.81
timestep = 0.02
friction_coefficient = 0.8
damper = 0.5
spring = 1.0
env = get_environment(:quadruped,
    representation=:minimal,
    timestep=timestep,
    contact_body=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring,
    infeasible_control=true,
    vis=vis)

# ## dimensions
nx = env.num_states
nu = env.num_inputs
nu_infeasible = 6

################################################################################
# Load trajectory and gains
################################################################################
file = JLD2.jldopen(joinpath(@__DIR__, "../data/trotting_forward.jld2"))
x_ref = file["x"]
u_ref = file["u"]
JLD2.close(file)

file = JLD2.jldopen(joinpath(@__DIR__, "../data/tvlqr_gains.jld2"))
K_tvlqr = file["K"]
JLD2.close(file)

file = JLD2.jldopen(joinpath(@__DIR__, "../data/tuned_gains.jld2"))
K_tuned = file["K"]
JLD2.close(file)

T = length(x_ref)


################################################################################
# Simulate accurately the TVLQR policy
################################################################################
# divide timestep
S = 1
mech_sim = get_mechanism(:quadruped,
    contact_body=true,
    timestep=timestep/S,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)

# controller
function ctrl!(m, k; x_sol=x_sol, u_sol=u_sol, K_sol=K_sol)
    N = length(u_sol)
    nu = input_dimension(m)
    @show Int(floor((k-1)/S))
    ind = (Int(floor((k-1)/S))) % N + 1

    x = get_minimal_state(m)
    # x[1] -= 0.15
    u = u_sol[ind] - K_sol[ind] * (x - x_sol[ind])
    u = clamp.(u, -0.2, +0.2)
    u = [zeros(6); u[7:end]] / S
    set_input!(m, SVector{nu}(u))
end

controller_tvlqr!(m, k) = ctrl!(m, k; K_sol=K_tvlqr)
controller_tuned!(m, k) = ctrl!(m, k; K_sol=K_tuned)

x1_roll = deepcopy(xref[1])
x1_roll[3] += 0.25
x1_roll[4] += 0.0
set_minimal_state!(mech_sim, x1_roll)
Main.@elapsed storage_tvlqr = simulate!(mech_sim, 4.0, controller_tvlqr!,
    record=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0),
    )

x1_roll = deepcopy(xref[1])
x1_roll[3] += 0.25
x1_roll[4] += 0.0
set_minimal_state!(mech_sim, x1_roll)
Main.@elapsed storage_tuned = simulate!(mech_sim, 4.0, controller_tuned!,
    record=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0),
    )

vis, anim = Dojo.visualize(mech_sim, storage_tvlqr, vis=vis, build=false, name=:tvlqr)
vis, anim = Dojo.visualize(mech_sim, storage_tuned, vis=vis, animation=anim, build=false, name=:tuned)

settransform!(vis[:tuned], Translation(0,0.5,0))
# convert_frames_to_video_and_gif("quadruped_tuned_vs_tvlqr_facing")
