# PREAMBLE

# PKG_SETUP

# ## Setup
using Dojo 

# ## Visualizer
vis = Visualizer()
open(vis)

################################################################################
## Astronaut
################################################################################
## humanoid
## multiple bodies
## initial linear velocities
## no gravity
## no spring and damper
## no control
################################################################################

function astronaut_simulation(vis::Visualizer; 
        timestep=1.0e-2, 
        tsim=1.0, 
        gravity=0.0,
        spring=0.0, 
        damper=0.0, 
        seed=100, 
        ϵ=1.0e-14, 
        display=true)

    Random.seed!(seed)
    mech = get_mechanism(:humanoid, 
        timestep=timestep, 
        gravity=gravity, 
        spring=spring, 
        damper=damper, 
        contact=false)
    initialize!(mech, :humanoid)

    ## Initialize bodies with random velocities
    for body in mech.bodies
        v = rand(3) .- 0.5
        v ./= norm(v)
        set_maximal_velocities!(body, 
            v=v)
    end

    tcompute = @elapsed storage = simulate!(mech, tsim, 
        record=true, 
        opts=SolverOptions(rtol=ϵ, btol=ϵ))

    return mech, storage, tcompute
end

# ## Simulation
solver_tolerance = 1.0e-14
timestep = [0.05, 0.01, 0.005]
tsim = 10.0
storage = []
mech, _ = astronaut_simulation(vis, 
    tsim=tsim, 
    timestep=timestep[1], 
    ϵ=solver_tolerance)

for t in timestep
    mech, data, tcompute = astronaut_simulation(vis, 
        tsim=tsim, 
        timestep=t, 
        ϵ=solver_tolerance)
    push!(storage, data)
end


# ## Energy
color = [:magenta, :orange, :cyan];
me = [Dojo.mechanical_energy(mech, s)[2:end] for s in storage]
plt = plot(xlabel="time (s)", ylabel="energy drift");
for (i, e) in enumerate(me)
    tt = [j * timestep[i] for j = 1:length(e)]
    plt = plot!(tt, e .- e[1], color=color[i], width=2.0, label="h = $(timestep[i])")
end
display(plt);

# ## Momentum
m = [Dojo.momentum(mech, s)[1:end] for s in storage]

## linear
mlin = [[Vector(mi)[1:3] for mi in mt] for mt in m]
plt = plot(ylims=(-1.0e-10, 1.0e-10), xlabel="time (s)", ylabel="linear momentum drift");
for (i, mt) in enumerate(mlin)
    tt = [j * timestep[i] for j in 1:length(mt)]
    plt = plot!(tt, hcat([mi - mt[1] for mi in mt]...)', color=color[i], label=["h = $(timestep[i])" "" ""])
end
display(plt)

## angular 
mang = [[Vector(mi)[4:6] for mi in mt] for mt in m]
plt = plot(ylims=(-1.0e-10, 1.0e-10), xlabel="time (s)", ylabel="angular momentum drift");
for (i, mt) in enumerate(mang)
    tt = [j * timestep[i] for j in 1:length(mt)]
    plt = plot!(tt, hcat([mi - mt[1] for mi in mt]...)', color=color[i], label=["h = $(timestep[i])" "" ""])
end
display(plt)