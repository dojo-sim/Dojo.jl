using Polyhedra

vis = Visualizer()
open(vis)

# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

include("methods.jl")

file = jldopen(joinpath(module_dir(), "examples",
	"real2sim", "hardware_examples", "sol_best6.jld2"))
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



vis = Visualizer()
open(vis)

vis, anim = cone_morphing(Dsol, vis=vis, fps=20, rot=0.03,
	vis_truth=true, vis_learned=true, translate=true, cam_pos=[0,-2,0.7], alt=-1.5)

t = 1
cam_pos = [0,-1.5,0.5]
vis, anim = cone_morphing(Dsol[t:t], vis=vis, fps=20, rot=0.00, background=true,
	vis_truth=false, vis_learned=true, translate=false, cam_pos=cam_pos, alt=-1.65, b0=0, b1=0)
