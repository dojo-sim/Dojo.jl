using PGFPlots
const PGF = PGFPlots

color = ["magenta", "orange", "cyan"];
p_energy = [PGF.Plots.Linear([j * timestep[i] for j = 1:length(e)], e .- e[1], legendentry="h=$(timestep[i])", mark="none",style="color="*color[i]*", line width=2pt, solid") for (i, e) in enumerate(me)]

a1 = Axis(vcat(p_energy...),
    xmin=0,
    xmax=10.0,
    hideAxis=false,
	ylabel="energy drift",
	xlabel="time (s)",
    legendPos="south east")
dir = joinpath("/home/taylor/Downloads")
PGF.save(joinpath(dir, "energy_drift.tikz"), a1)

##########
# Visualization
color = RGBA(255.0/255.0,0.0,255.0,1.0);
z = get_maximal_state(storage)
z = [[z[1] for t = 1:100]..., z..., [z[end] for t = 1:100]...]
build_robot(mech, vis=vis, color=color)
T = length(z)
anim = MeshCat.Animation(convert(Int, floor(1.0 / timestep0)))
for t = 1:T
    MeshCat.atframe(anim, t) do
        set_robot(vis, mech, z[t])
    end
end
MeshCat.setanimation!(vis, anim)
setvisible!(vis["/Axes"], false)
setvisible!(vis["/Grid"], false)

# Ghost
t = 400 # 225, 150, 1
set_robot(vis, mech, z[t])
##########