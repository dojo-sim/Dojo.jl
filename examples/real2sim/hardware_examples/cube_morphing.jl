using Polyhedra

vis = Visualizer()
open(vis)

# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

function cube_morphing(D; vis=Visualizer(), fps=30, rot=0.00, vis_truth::Bool=true,
		vis_learned::Bool=true, translate::Bool=true, cam_pos=[0,-10,4.], zoom=1.8,
		color_truth=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0),
		color_learn=RGBA(51.0 / 255.0, 1.0, 1.0, 1.0),
		background=true, alt=0.0, b0=0.2, b1=1.5)
	setvisible!(vis["/Background"], background)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	setvisible!(vis["/Grid"], false)
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", zoom)

	settransform!(vis["/Cameras/default/rotated/<object>"], MeshCat.Translation(0.,0, 0))
	settransform!(vis["/Cameras/default/rotated"], MeshCat.Translation(0.,0, 0))
	settransform!(vis["/Cameras/default"], MeshCat.Translation(0.,0, 0))
	settransform!(vis["/Cameras"], MeshCat.Translation(cam_pos...))

	anim = MeshCat.Animation(fps)
	b0 = Int(floor(b0*fps))
	b1 = Int(floor(b1*fps))
	D = [fill(D[1], b0); D; fill(D[end] .+ 0.005*rand(length(D[end])), 2b1)]
	setobject!(vis[:cube],
		Rect(Vec(-1,-1,-1),Vec(2,2,2.0)),
		MeshPhongMaterial(color=color_truth))

	for (i,d) in enumerate(D)
		v = vrep([d[1 + (i-1)*3 .+ (1:3)] for i = 1:8])
		p = polyhedron(v)
		Polyhedra.Mesh(p)
		setobject!(vis[Symbol("polyhedra$i")],
			Polyhedra.Mesh(p),
			MeshPhongMaterial(color=color_learn))
	end

	for (i,d) in enumerate(D)
		atframe(anim, i) do
			MeshCat.setvisible!(vis[:cube], vis_truth)
			for (j,d) in enumerate(D)
				MeshCat.setvisible!(vis[Symbol("polyhedra$j")], i==j && vis_learned)
			end
			settransform!(vis[:cube], MeshCat.compose(
				MeshCat.Translation(-2.0*translate,0,1+alt),
				MeshCat.LinearMap(RotZ(rot*min(i,length(D)-b1))),
				))

			settransform!(vis[Symbol("polyhedra$i")], MeshCat.compose(
				MeshCat.Translation(2.0*translate,0,1+alt),
				MeshCat.LinearMap(RotZ(rot*min(i,length(D)-b1))),
				))

		end
	end
	setanimation!(vis, anim)
	return vis, anim
end


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

function cube_sim_v_truth(d, traj_truth::Storage{T,HT}, traj_sim::Storage{T,HS};
		vis=Visualizer(), fps=30, zoom=1.8,
		transparency_truth=1.0, 
		transparency_learn=1.0,
		color_truth=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, transparency_truth),
		color_learn=RGBA(51.0 / 255.0, 1.0, 1.0, transparency_learn),
		background=true, b0=0.2, b1=1.5) where {T,HT,HS}
	setvisible!(vis["/Background"], background)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	# setvisible!(vis["/Grid"], false)
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", zoom)
	setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 0.80)


	anim = MeshCat.Animation(fps)
	b0 = Int(floor(b0*fps))
	b1 = Int(floor(b1*fps))
	H = max(HT,HS)

	# Create Objects
	setobject!(vis[:truth],
		Rect(Vec(-1,-1,-1),Vec(2,2,2.0)),
		MeshPhongMaterial(color=color_truth))

	v = vrep([d[1 + (i-1)*3 .+ (1:3)] for i = 1:8])
	p = polyhedron(v)
	Polyhedra.Mesh(p)
	setobject!(vis[:sim],
		Polyhedra.Mesh(p),
		MeshPhongMaterial(color=color_learn))

	# Animate
	for i = 1:b0
		atframe(anim, i) do
			settransform!(vis[:truth], MeshCat.compose(
				MeshCat.Translation(traj_truth.x[1][1]),
				MeshCat.LinearMap(traj_truth.q[1][1]),
				))
			settransform!(vis[:sim], MeshCat.compose(
				MeshCat.Translation(traj_sim.x[1][1]),
				MeshCat.LinearMap(traj_sim.q[1][1]),
				))
		end
	end
	for i = 1:H
		atframe(anim, b0+i) do
			settransform!(vis[:truth], MeshCat.compose(
				MeshCat.Translation(traj_truth.x[1][min(i,HT)]),
				MeshCat.LinearMap(traj_truth.q[1][min(i,HT)]),
				))
			settransform!(vis[:sim], MeshCat.compose(
				MeshCat.Translation(traj_sim.x[1][min(i,HS)]),
				MeshCat.LinearMap(traj_sim.q[1][min(i,HS)]),
				))
		end
	end
	for i = 1:b1
		atframe(anim, b0+H+i) do
			settransform!(vis[:truth], MeshCat.compose(
				MeshCat.Translation(traj_truth.x[1][end]),
				MeshCat.LinearMap(traj_truth.q[1][end]),
				))
			settransform!(vis[:sim], MeshCat.compose(
				MeshCat.Translation(traj_sim.x[1][end]),
				MeshCat.LinearMap(traj_sim.q[1][end]),
				))
		end
	end
	setanimation!(vis, anim)
	return vis, anim
end

function cube_ghost_sim_v_truth(d, traj_truth::Storage{T,HT}, traj_sim::Storage{T,HS};
		vis=Visualizer(), zoom=1.8,	
		transparency_truth=1.0, 
		transparency_learn=1.0,
		color_truth=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, transparency_truth),
		color_learn=RGBA(51.0 / 255.0, 1.0, 1.0, transparency_learn),
		background=true) where {T,HT,HS}
	setvisible!(vis["/Background"], background)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	# setvisible!(vis["/Grid"], false)
	setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", zoom)
	setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", 0.80)

	H = max(HT,HS)
	frames = [5, 30, H]
	α = [0.3, 0.5,1.0]
	# Create Objects
	for (j,i) in enumerate(frames)
		setobject!(vis["truth$(i)"],
			Rect(Vec(-1,-1,-1),Vec(2,2,2.0)),
			MeshPhongMaterial(color=color_truth))

		v = vrep([d[1 + (i-1)*3 .+ (1:3)] for i = 1:8])
		p = polyhedron(v)
		Polyhedra.Mesh(p)
		setobject!(vis["sim$(i)"],
			Polyhedra.Mesh(p),
			MeshPhongMaterial(color=color_learn))
	end

	# Animate
	for i in frames
		settransform!(vis["truth$(i)"], MeshCat.compose(
			MeshCat.Translation(traj_truth.x[1][min(i,HT)]),
			MeshCat.LinearMap(traj_truth.q[1][min(i,HT)]),
			))
		settransform!(vis["sim$(i)"], MeshCat.compose(
			MeshCat.Translation(traj_sim.x[1][min(i,HS)]),
			MeshCat.LinearMap(traj_sim.q[1][min(i,HS)]),
			))
	end
	return vis
end
