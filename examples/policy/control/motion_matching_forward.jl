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
using Plots

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
# Load trajectory
################################################################################
velocities = [0, 0.02, 0.04, 0.07, 0.10, 0.13, 0.16, 0.19, 0.22, 0.25]
x_lib = Dict()
u_lib = Dict()

for velocity in velocities
    file = JLD2.jldopen(joinpath(@__DIR__, "../data_forward/trotting_variable_$(velocity).jld2"))
    x_lib[velocity] = file["x"]
    u_lib[velocity] = file["u"]
    JLD2.close(file)
end
T = length(x_lib[0])

################################################################################
# TVLQR
################################################################################
K_lib = Dict()
q_tracking=[0.3; 0.05; 0.05;
    5e-1 * ones(3);
    1e-6 * ones(3);
    1e-6 * ones(3);
    fill([2, 1e-6], 12)...]
r_tracking = env.mechanism.timestep * 30 * ones(nu)
for velocity in velocities
    K_lib[velocity] = tvlqr(x_lib[velocity], u_lib[velocity], env;
        q_tracking=q_tracking, r_tracking=r_tracking)[1]
end

################################################################################
# Target trajectory
################################################################################
# divide timestep
S = 3
# v_target = [0.20 * (1 - cos(2π * t/15))/2 for t=0:timestep/S:15]
v_target = [0.20 * (1 - cos(2π * t/7))/2 for t=0:timestep/S:15]
plot(0:timestep/S:15, v_target)

################################################################################
# Simulate accurately the TVLQR policy
################################################################################
mech_sim = get_mechanism(:quadruped,
    contact_body=true,
    timestep=timestep/S,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)


function find_closest(velocities, v)
    Δ = abs.(velocities .- v)
    velocity = velocities[findmin(Δ)[2]]
    return velocity
end

# controller
function ctrl!(mechanism, k; velocities=velocities, x_lib=x_lib, u_lib=u_lib, K_lib=K_lib)
    @show k
    # choose reference based on target velocity
    velocity = find_closest(velocities, v_target[k])
    x_ref = x_lib[velocity]
    u_ref = u_lib[velocity]
    K_ref = K_lib[velocity]

    # find index in the loops
    N = length(u_ref)
    nu = input_dimension(mechanism)
    ind = (Int(floor((k-1)/S))) % N + 1

    # compute policy
    x = get_minimal_state(mechanism)
    Δx = (x_ref[ind] - x)
    Δx[1] = 0.0
    u = u_ref[ind] + K_ref[ind] * Δx
    u = clamp.(u, -0.2, +0.2)
    u = [zeros(6); u[7:end]] / S
    set_input!(mechanism, SVector{nu}(u))
end

controller!(m, k) = ctrl!(m, k; velocities=velocities, x_lib=x_lib, u_lib=u_lib, K_lib=K_lib)

# x1_roll = deepcopy(x_lib[0.10][1])
x1_roll = deepcopy(x_lib[0.0][1])
x1_roll[3] += 0.0
x1_roll[4] += 0.0
set_minimal_state!(mech_sim, x1_roll)
Main.@elapsed storage = simulate!(mech_sim, 15.0, controller!,
    record=true,
    opts=SolverOptions(rtol=1e-5, btol=1e-4, undercut=5.0),
    )

vis, anim = Dojo.visualize(mech_sim, storage, vis=vis, build=true, name=:robot)
# convert_frames_to_video_and_gif("quadruped_motion_matching")
