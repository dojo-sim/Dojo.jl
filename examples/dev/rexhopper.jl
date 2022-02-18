using Dojo 

# Open visualizer
vis = Visualizer()
# open(vis)
render(vis)

# Mechanism
mechanism = get_rexhopper(model=:rexhopper2, timestep=0.01, gravity=-9.81, contact_body=true, friction_coefficient=1.0)

for body in mechanism.bodies 
    body.inertia = 10.0 * body.inertia
    @show body.inertia 
    @show norm(diag(body.inertia))
end

# Simulate
initialize!(mechanism, :rexhopper, x=[0.0; 0.0; 0.4])

# Open visualizer
storage = simulate!(mechanism, 5.0, record=true, opts=SolverOptions(undercut=10.0, btol=1.0e-4, rtol=1.0e-4, verbose=true));

# Visualize
visualize(mechanism, storage, vis=vis, show_contact=true);
