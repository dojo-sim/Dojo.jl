# Defining a Simulation
Here, we explain how to simulate a dynamical system i.e., a [`Mechanism`](@ref) forward in time.
The example that we are trying to replicate the Dzhanibekov effect shown below.

![dzhanibekov](./assets/dzhanibekov_nasa.gif)

Load the `Dojo` package.
```julia
using Dojo
```

Define the simulation time step, and other parameters like the gravity.
```julia
timestep=0.01
gravity=0.0
```

We use the mechanism called `:dzhanibekov`.
```julia
mech = get_mechanism(:dzhanibekov,
        timestep=timestep,
        gravity=gravity);
```

We initialize the system with a given initial angular velocity.
```julia
initialize!(mech, :dzhanibekov,
    angular_velocity=[15.0; 0.01; 0.0])
```

We simulate this system for 5 seconds, we record the resulting trajectory in `storage`,
```julia
storage = simulate!(mech, 5.0,
    record=true,
    verbose=false)
```

We visualize the trajectory in the browser,
```julia
vis = Visualizer()
open(vis)
visualize(mech, storage, vis=vis)
```

And voila! You should see something like this;

```@raw html
<img src="./assets/animations/dzhanibekov.gif" width="300"/>
```
