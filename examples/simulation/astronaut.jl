using Dojo 

# visualizer
vis = Visualizer()
open(vis)

################################################################################
# Astronaut
################################################################################
# humanoid
# multiple bodies
# initial linear velocities
# no gravity
# no spring and damper
# no control
################################################################################

function astronaut_simulation(vis::Visualizer; Δt=1e-2, tsim=1.0, g=0.0,
        spring=0.0, damper=0.0, seed::Int=100, ϵ=1e-14, display::Bool = true)
    
    Random.seed!(seed)
    mech = getmechanism(:humanoid, Δt=Δt, g=g, spring=spring, damper=damper, contact=false)
    initialize!(mech, :humanoid)

    # Initialize bodies with random velocities
    for body in mech.bodies
        v = rand(3) .- 0.5
        v ./= norm(v)
        setVelocity!(body, v = v)
    end

    tcompute = @elapsed storage = simulate!(mech, tsim, record=true, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))

    return mech, storage, tcompute
end

# Simulation
solver_tolerance = 1.0e-14
timestep = [0.05, 0.01, 0.005]
tsim = 10.0
storage = Storage[]
mech, _ = astronaut_simulation(vis, tsim=tsim, Δt=timestep[1], ϵ=solver_tolerance)
for t in timestep
    mech, data, tcompute = astronaut_simulation(vis, tsim=tsim, Δt=t, ϵ=solver_tolerance)
    push!(storage, data)
end

# visualize(mech, storage, vis=vis)
# color_magenta = RGBA(255.0 / 255.0, 0.0, 255.0, 1.0);
# color_orange = RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0);
# color_cyan = RGBA(51.0 / 255.0, 1.0, 1.0, 1.0);
color = [:magenta, :orange, :cyan];

# energy
me = [Dojo.mechanicalEnergy(mech, s)[2:end] for s in storage]
plt = plot(xlabel="time (s)", ylabel="energy drift");
for (i, e) in enumerate(me)
    tt = [j * timestep[i] for j = 1:length(e)]
    plt = plot!(tt, e .- e[1], color=color[i], width=2.0, label="h = $(timestep[i])")
end 
display(plt);

# momentum
m = [Dojo.momentum(mech, s)[1:end] for s in storage]

# linear momentum
mlin = [[Vector(mi)[1:3] for mi in mt] for mt in m]
plt = plot(ylims=(-1.0e-10, 1.0e-10), xlabel="time (s)", ylabel="linear momentum drift");
for (i, mt) in enumerate(mlin)
    tt = [j * timestep[i] for j in 1:length(mt)]
    plt = plot!(tt, hcat([mi - mt[1] for mi in mt]...)', color=color[i], label=["h = $(timestep[i])" "" ""])
end
display(plt)

# angular momentum
mang = [[Vector(mi)[4:6] for mi in mt] for mt in m]
plt = plot(ylims=(-1.0e-10, 1.0e-10), xlabel="time (s)", ylabel="angular momentum drift");
for (i, mt) in enumerate(mang)
    tt = [j * timestep[i] for j in 1:length(mt)]
    plt = plot!(tt, hcat([mi - mt[1] for mi in mt]...)', color=color[i], label=["h = $(timestep[i])" "" ""])
end
display(plt)

# 
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
z = getMaxState(storage)
z = [[z[1] for t = 1:100]..., z..., [z[end] for t = 1:100]...]
build_robot(vis, mech, color=color)
T = length(z) 
anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt0)))
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

