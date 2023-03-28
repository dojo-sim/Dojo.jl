# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
# Pkg.instantiate()

# ## Setup
using Dojo

# ## Mechanism
mech = get_mechanism(:tippetop,
    timestep=0.01,
    gravity=-9.81,
    friction_coefficient=0.4,
    contact=true,
    contact_type=:nonlinear)

# ## Simulate
initialize!(mech, :tippetop,
    body_position=[0.0, 0.0, 1.0],
    body_orientation=[0,0.01,0],
    body_angular_velocity=[0.0, 0.01, 50.0])

storage = simulate!(mech, 10.0,
    record=true,
    verbose=false,
    opts=SolverOptions(verbose=false, btol=1e-5))

# ## Open visualizer
vis = Visualizer()
render(vis)
visualize(mech, storage,
    vis=vis)
