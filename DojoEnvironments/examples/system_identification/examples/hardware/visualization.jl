using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
Pkg.instantiate()

# ## Setup
using Dojo
using Plots
using Random
using MeshCat
using Polyhedra
using JLD2
using DojoEnvironments

# ## Include methods
methods_dir = joinpath("../../methods")
include(joinpath(methods_dir, "hardware_visualization.jl"))


# ## ---------------------------------------------------------------------------
# parameters
# ## ---------------------------------------------------------------------------
S = 1
# we multiplied all the distance by length_scaling for better data scaling
# we must rescale the gravity term accordingly since gravity contains length (g == kg.m.s^-2)
length_scaling = 20.0
timestep = 1/148 * S
gravity = -9.81 * scaling
friction_coefficient = 0.16
radius = 0.00
side = 2.00
mode = :box
opts_step = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
opts_grad = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
model = :block
N = 100

mech_kwargs = Dict(
	:timestep => timestep,
	:gravity => gravity,
	:friction_coefficient => friction_coefficient,
	:radius => radius,
	:side => side,
	:mode => mode)

# ## ---------------------------------------------------------------------------
# Load Dataset
# ## ---------------------------------------------------------------------------
params0, trajs0 = open_dataset(model; experiment_type="hardware", N=N, mech_kwargs...)
data0 = params0[:data]
nd = sum(data_dim.(mech.contacts))
data_contacts0 = data0[end-nd+1:end]

# ## ---------------------------------------------------------------------------
# Load solution
# ## ---------------------------------------------------------------------------
file = jldopen(joinpath(@__DIR__, "../..", "data", "hardware", "visualization",
	"hardware_block_solution.jld2"))
dsol = file["dsol"]


# ## ---------------------------------------------------------------------------
# We learn a single coefficient of friction and a single side length -> 2 params in total
# ## ---------------------------------------------------------------------------
function d_to_data_contacts(d)
	friction_coefficient = d[1]
	half_side = d[2]
	data_contacts = [
			friction_coefficient; 0.0; +half_side; +half_side; -half_side;
			friction_coefficient; 0.0; +half_side; -half_side; -half_side;
			friction_coefficient; 0.0; -half_side; +half_side; -half_side;
			friction_coefficient; 0.0; -half_side; -half_side; -half_side;
			friction_coefficient; 0.0; +half_side; +half_side; +half_side;
			friction_coefficient; 0.0; +half_side; -half_side; +half_side;
			friction_coefficient; 0.0; -half_side; +half_side; +half_side;
			friction_coefficient; 0.0; -half_side; -half_side; +half_side;
			]
	return data_contacts
end

data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(2))


mech = get_mechanism(model; mech_kwargs...)
set_data!(mech.contacts, d_to_data_contacts(dsol[1]))

################################################################################
# # Learned vs truth: trajectory
################################################################################

# ## Open visualizer
vis = Visualizer()
render(vis)

id = 7
traj_data = trajs0[id]
z0 = get_maximal_state(traj_data, 1)
set_maximal_state!(mech, z0)

traj_sim = simulate!(mech, 0.80, opts=opts_step)

cube_sim_v_truth(dsol[1], traj_data, traj_sim, vis=vis,
	transparency_truth=0.5,
	fps=Int(floor(1/mech.timestep)), b0=0.0, b1=0.0)

cube_sim_v_truth(dsol[1], traj_data, traj_sim, vis=vis, transparency_truth=1.0)
cube_ghost_sim_v_truth(dsol[1], traj_data, traj_sim, vis=vis, transparency_truth=1.0)
set_floor!(vis, x=20, y=20, color=RGBA(ones(4)...))
set_light!(vis, ambient=0.80)

set_camera!(vis, zoom=4.0)


################################################################################
# # Learned vs truth: geometry and friction cone
################################################################################

# ## Open visualizer
vis = Visualizer()
render(vis)

cam_pos = [2,-4.5,1.8]
vis, anim = cube_morphing(dsol[2], vis=vis, fps=20, rot=0.03,
	vis_truth=true, vis_learned=true, translate=true, cam_pos=[0,-6,1.8], alt=-1)
vis, anim = cube_morphing(dsol[2][5:5], vis=vis, fps=20, rot=0.00, background=false,
	vis_truth=false, vis_learned=true, translate=false, cam_pos=cam_pos, alt=-1, b0=0, b1=0)
vis, anim = cube_morphing(dsol[2], vis=vis, fps=20, rot=0.00, background=false,
	vis_truth=true, vis_learned=false, translate=false, cam_pos=cam_pos, alt=-1, b0=0, b1=0)

vis = Visualizer()
render(vis)
t = 21
vis, anim = cube_morphing(dsol[2][t:t], vis=vis, fps=20, rot=0.00, background=false,
	vis_truth=true, vis_learned=true, translate=false, cam_pos=cam_pos, alt=-1, b0=0, b1=0)

vis = Visualizer()
render(vis)

vis, anim = cone_morphing(dsol[2], vis=vis, fps=20, rot=0.03,
	vis_truth=true, vis_learned=true, translate=true, cam_pos=[0,-15,2.0], alt=-1.5, zoom=20.0)

t = 1
cam_pos = [0,-1.5,0.5]
vis, anim = cone_morphing(dsol[2][t:t], vis=vis, fps=20, rot=0.00, background=true,
	vis_truth=false, vis_learned=true, translate=false, cam_pos=cam_pos, alt=-1.65, b0=0, b1=0)
