# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

# ## Setup
using Dojo 
using DojoEnvironments
using Random
using LinearAlgebra
using Plots

seed=100
Random.seed!(seed)

function simulation_loop(timestep, tsim;
        gravity=0.0, springs=0.0, dampers=0.0, ϵ=1.0e-6)

    display(timestep)

    mech = get_mechanism(:humanoid; 
        timestep, gravity, springs, dampers,
        contact_feet=false)

    ## Initialize bodies with random velocities
    for body in mech.bodies
        v = LinearAlgebra.normalize(randn(3))
        set_maximal_velocities!(body; v=v)
    end

    storage = simulate!(mech, tsim; 
        record=true, opts=SolverOptions(btol=ϵ))

    return mech, storage
end

# ## Simulation
timesteps = [0.05, 0.01, 0.005]
tsim = 10.0
storages = []

for timestep in timesteps
    mech, storage = simulation_loop(timestep, tsim)
    visualize(mech, storage, visualize_floor=false)
    push!(storages, storage)
end

# ## Plots
color = [:magenta, :orange, :cyan];
mech, _ = simulation_loop(0.01, 0.01)
conservation_plot = plot(layout=(3,1), size=(800,500), legend=:topleft)

# ## Energy
energies = [Dojo.mechanical_energy(mech, storage)[2:end] for storage in storages]
energy_plot = plot(xlabel="time (s)", ylabel="energy drift");
for (i, energy) in enumerate(energies)
    time = [step * timesteps[i] for step = 1:length(energy)]
    plot!(conservation_plot[1,1], time, energy .- energy[1]; width=2.0, color=color[i],
        xlabel="time (s)", ylabel="energy drift", label="h = $(timesteps[i])")
    display(conservation_plot)
end

# ## Momentum
momenta = [Dojo.momentum(mech, storage)[1:end] for storage in storages]

# ## linear
momenta_linear = [[Vector(momentum_i)[1:3] for momentum_i in momentum] for momentum in momenta]
for (i, momentum) in enumerate(momenta_linear)
    time = [step * timesteps[i] for step = 1:length(momentum)]
    plot!(conservation_plot[2,1], time, hcat([momentum_i - momentum[1] for momentum_i in momentum]...)'; width=2.0, color=color[i],
        xlabel="time (s)", ylabel="linear mom. drift", label=["h = $(timesteps[i])" "" ""])
    display(conservation_plot)
end

# ## angular 
momenta_angular = [[Vector(momentum_i)[4:6] for momentum_i in momentum] for momentum in momenta]
momentum_angular_plot = plot(ylims=(-1.0e-10, 1.0e-10), xlabel="time (s)", ylabel="angular momentum drift");
for (i, momentum) in enumerate(momenta_angular)
    time = [j * timesteps[i] for j in 1:length(momentum)]
    plot!(conservation_plot[3,1], time, hcat([momentum_i - momentum[1] for momentum_i in momentum]...)'; width=2.0, color=color[i],
        xlabel="time (s)", ylabel="angular mom. drift", label=["h = $(timesteps[i])" "" ""])
    display(conservation_plot)
end