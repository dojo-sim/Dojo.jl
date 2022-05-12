using Dojo
using Polyhedra
using JLD2


################################################################################
# # Learned vs truth: geometry and friction cone
################################################################################
include("methods.jl")

vis= Visualizer()
open(vis)

file = jldopen(joinpath(@__DIR__, "..", "results", "sol_best6.jld2"))
Dsol = file["Dsol"]
cam_pos = [2,-4.5,1.8]
vis, anim = cube_morphing(Dsol, vis=vis, fps=20, rot=0.03,
	vis_truth=true, vis_learned=true, translate=true, cam_pos=[0,-6,1.8], alt=-1)
vis, anim = cube_morphing(Dsol[5:5], vis=vis, fps=20, rot=0.00, background=false,
	vis_truth=false, vis_learned=true, translate=false, cam_pos=cam_pos, alt=-1, b0=0, b1=0)
vis, anim = cube_morphing(Dsol, vis=vis, fps=20, rot=0.00, background=false,
	vis_truth=true, vis_learned=false, translate=false, cam_pos=cam_pos, alt=-1, b0=0, b1=0)

t = 50
vis, anim = cube_morphing(Dsol[t:t], vis=vis, fps=20, rot=0.00, background=false,
	vis_truth=true, vis_learned=false, translate=false, cam_pos=cam_pos, alt=-1, b0=0, b1=0)

vis= Visualizer()
open(vis)

vis, anim = cone_morphing(Dsol, vis=vis, fps=20, rot=0.03,
	vis_truth=true, vis_learned=true, translate=true, cam_pos=[0,-15,2.0], alt=-1.5, zoom=20.0)

t = 1
cam_pos = [0,-1.5,0.5]
vis, anim = cone_morphing(Dsol[t:t], vis=vis, fps=20, rot=0.00, background=true,
	vis_truth=false, vis_learned=true, translate=false, cam_pos=cam_pos, alt=-1.65, b0=0, b1=0)


################################################################################
# # Learned vs truth: trajectory
################################################################################
include(joinpath(@__DIR__, "..", "methods", "utils.jl"))
include(joinpath(@__DIR__, "methods.jl"))

# open visualizer
vis= Visualizer()
open(vis)

S = 7
timestep= 1/148 * S
gscaled = -9.81*20

file = jldopen(joinpath(@__DIR__, "..", "results", "sol_best6.jld2"))
Dsol = file["Dsol"]

# Load Dataset
params0, trajs0, pairs0 = open_dataset(:hardwarebox; N=400, s=S)
params1, trajs1, pairs1 = open_dataset(:hardwarebox; N=400, s=1)


function d2data(d)
	friction_coefficient = d[1]
	data = [friction_coefficient; 0;0;0; +d[2:4];
			friction_coefficient; 0;0;0; +d[5:7];
			friction_coefficient; 0;0;0; +d[8:10];
			friction_coefficient; 0;0;0; +d[11:13];
			friction_coefficient; 0;0;0; +d[14:16];
			friction_coefficient; 0;0;0; +d[17:19];
			friction_coefficient; 0;0;0; +d[20:22];
			friction_coefficient; 0;0;0; +d[23:25];
			]
	return data
end
mech = get_mechanism(:block, timestep=timestep/S, gravity=gravityscaled, friction_coefficient=Dsol[end][1], radius=0.00, side=2.0, mode=:box);
set_simulator_data!(mech, d2data(Dsol[end]))
id = 7#4,6,7,8
traj_truth = trajs1[id]
x2 = traj_truth.x[1][1] - [0,0,2.0]/2
v15 = traj_truth.v[1][1]
q2 = traj_truth.q[1][1]
ϕ15 = traj_truth.ω[1][1]

initialize!(mech, :block, x=x2, v=v15, q=q2, ω=ϕ15)
traj_sim = simulate!(mech, 0.80, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))

cube_sim_v_truth(Dsol[end], traj_truth, traj_sim, vis=vis,
	transparency_truth=0.5,
	fps=Int(floor(1/mech.timestep)), b0=0.0, b1=0.0)

cube_sim_v_truth(Dsol[end], traj_truth, traj_sim, vis=vis, transparency_truth=1.0)
cube_ghost_sim_v_truth(Dsol[end], traj_truth, traj_sim, vis=vis, transparency_truth=1.0)
set_floor!(vis, x=20, y=20, color=RGBA(ones(4)...))
set_light!(vis, ambient=0.80)

set_camera!(vis, zoom=4.0)
