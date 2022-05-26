using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo

# ## Mechanism
mech = get_mechanism(:hopper, 
    timestep=0.1, 
    gravity=-9.81, 
    # friction_coefficient=0.5, 
    damper=1.0, 
    spring=0.0, 
    contact_foot=true, 
    contact_body=true, 
   )

# ## Simulate
Dojo.initialize!(mech, :hopper)

storage = simulate!(mech, 2.5, 
    record=true, 
    opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6, verbose=true))

# ## Visualize
vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis, show_contact=true)

# ## Contact interpenetration 
res = get_sdf(mech, storage) # distance from floor to each contact
minimum(minimum([min.(0.0, r) for r in res]))

