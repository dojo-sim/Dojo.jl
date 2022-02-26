# PREAMBLE

# PKG_SETUP

# ## Setup
using Dojo

# ## Simulation
timestep = 0.01
gravity = 0.0
mech = get_dzhanibekov(
        timestep=timestep, 
        gravity=gravity);

# ## Simulate
initialize_dzhanibekov!(mech, 
    ω=[15.0; 0.01; 0.0])
storage = simulate!(mech, 4.65, 
    record=true, verbose=false)

# ## Visualizers
vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis)

