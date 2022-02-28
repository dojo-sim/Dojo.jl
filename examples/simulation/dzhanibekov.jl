# PREAMBLE

# PKG_SETUP

# ## Setup
using Dojo

# ## Simulation
timestep=0.01
gravity=0.0
mech = get_mechanism(:dzhanibekov,
        timestep=timestep,
        gravity=gravity);

# ## Simulate
initialize_dzhanibekov!(mech,
    angular_velocity=[8.0; 0.01; 0.0]);
storage = simulate!(mech, 8.00,
    record=true,
    verbose=false);

# ## Visualizers
vis = Visualizer()
render(vis)
visualize(mech, storage, vis=vis);
