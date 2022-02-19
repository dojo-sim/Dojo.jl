using Dojo 

# Open visualizer
vis = Visualizer()
# open(vis)
render(vis)

# Mechanism
mechanism = get_rexhopper(model=:rexhopper_fixed, timestep=0.05, gravity=-9.81, contact_body=true, friction_coefficient=1.0)
env = make("rexhopper",
    timestep=0.05, gravity=-9.81, contact_body=true, friction_coefficient=1.0)

# Simulate
initialize!(env.mechanism, :rexhopper, x=[0.0; 0.0; 0.4])

# Open visualizer
storage = simulate!(mechanism, 2.5, record=true, opts=SolverOptions(undercut=10.0, btol=1.0e-4, rtol=1.0e-4, verbose=false));

# Visualize
visualize(mechanism, storage, vis=vis, show_contact=false);
