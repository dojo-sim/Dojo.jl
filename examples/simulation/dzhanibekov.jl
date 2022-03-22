using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo

# ## Simulation
timestep = 0.01
gravity = 0.0
mech = get_mechanism(:dzhanibekov,
        timestep=timestep,
        gravity=gravity);

# ## Simulate
initialize!(mech, :dzhanibekov,
    angular_velocity=[8.0; 0.01; 0.0]);
storage = simulate!(mech, 8.0,
    record=true,
    verbose=true)

# ## Visualizers
vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis);
