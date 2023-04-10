# ## Setup
using Dojo
using DojoEnvironments

# ## Get mechanism (check DojoEnvironment/mechanisms files for kwargs)
mechanism = get_mechanism(:pendulum; timestep=0.02, length=0.75) 

# ## Initialize mechanism (check DojoEnvironment/mechanisms files for kwargs)
initialize!(mech, :pendulum; angle=-0.3)

# ## Simulate mechanism
storage = simulate!(mechanism, 5, record=true)
    
# ## Visualize mechanism
visualize(mechanism, storage)
render(vis)