# PREAMBLE

# PKG_SETUP

# ## Setup
using Dojo

# ## Mechanism
mech = get_mechanism(:atlas, 
    timestep=0.1, 
    gravity=-9.81, 
    friction_coefficient=0.5, 
    damper=25.0, 
    spring=1.0, 
    contact_feet=true, 
    contact_body=true, 
    model_type=:v2)

# ## Simulate
initialize_atlasstance!(mech, 
    body_position=[0.0, 0.0, 0.5], 
    body_orientation=[0.0, 0.2, 0.1])

storage = simulate!(mech, 5.0, 
    record=true, 
    opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6, verbose=true))

# ## Visualize
vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis, show_contact=true)

# ## Contact innerpenetration 
res = get_sdf(mech, storage)
minimum(minimum([min.(0.0, r) for r in res]))

