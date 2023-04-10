# ## Setup
using Dojo
using DojoEnvironments

# ## Get mechanism (check DojoEnvironment/environments files for kwargs)
environment = get_environment(:pendulum; timestep=0.01, horizon=200) 

# ## Initialize mechanism (check DojoEnvironment/mechanisms files for kwargs)
initialize!(environment.mechanism, :pendulum; angle=-0.3)

# ## Simulate mechanism
simulate!(environment, record=true)
    
# ## Visualize mechanism
visualize(environment)
render(vis)