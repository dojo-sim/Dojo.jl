# ### Setup
# PKG_SETUP
using Dojo
using DojoEnvironments

# ### Get environment (check DojoEnvironment/environments files for kwargs)
environment = get_environment(:pendulum; timestep=0.01, horizon=200) 

# ### Initialize mechanism (check DojoEnvironment/mechanisms files for kwargs)
initialize!(environment.mechanism, :pendulum; angle=-0.3)

# ### Simulate environment
simulate!(environment, record=true)
    
# ### Visualize environment
vis = visualize(environment)
render(vis)