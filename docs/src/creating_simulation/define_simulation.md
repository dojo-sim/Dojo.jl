# Defining a Simulation
Here, we explain how to simulate a dynamical system i.e., a [`Mechanism`](@ref) forward in time.
The example that we are trying to replicate the Dzhanibekov effect shown below.

![dzhanibekov](../assets/dzhanibekov_nasa.gif)

```julia
# ### Setup
using Dojo
using DojoEnvironments

# ### Get mechanism (check DojoEnvironment/mechanisms files for kwargs)
mechanism = get_mechanism(:dzhanibekov; timestep=0.01, gravity=0) 

# ### Initialize mechanism (check DojoEnvironment/mechanisms files for kwargs)
initialize!(mechanism, :dzhanibekov; angular_velocity=[15.0; 0.01; 0.0])

# ### Simulate mechanism
storage = simulate!(mechanism, 5, record=true)
    
# ### Visualize mechanism
vis = visualize(mechanism, storage)
render(vis)
```

And voila! You should see something like this;

```@raw html
<img src="../assets/animations/dzhanibekov.gif" width="300"/>
```
