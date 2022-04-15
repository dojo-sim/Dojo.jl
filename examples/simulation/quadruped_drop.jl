using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo

# ## Mechanism
gravity = -9.81
timestep = 0.05
friction_coefficient = 0.8
damper = 1.0
env = get_environment(:quadruped,
    representation=:minimal,
    timestep=timestep,
    contact_body=true,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper);

# ## Simulate
initialize!(env.mechanism, :quadruped, 
    body_position=[0.0; 0.0; 0.5])

storage = simulate!(env.mechanism, 2.0, 
    record=true, 
    opts=SolverOptions(rtol=1.0e-4, btol=1.0e-4, verbose=false));

# ## Visualizer
open(env.vis)

# ## Visualize
visualize(env.mechanism, storage, 
    vis=env.vis, 
    show_contact=true)
