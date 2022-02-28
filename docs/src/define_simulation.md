# Defining a Simulation
Here, we explain how to simulate a dynamical system i.e. a [`mechanism`](@ref) forward in time.
The example that we are trying to replicate the Dzhanibekov effect shown below.

![dzhanibekov](./assets/dzhanibekov_nasa.gif)


Load the `Dojo` package.
```julia
using Dojo
```

Define the simulation timestep, and other parameters like the gravity.
```julia
timestep=0.01
gravity=0.0
```

We want to simulate the
```julia
mech = get_mechanism(:dzhanibekov,
        timestep=timestep,
        gravity=gravity);

# ## Simulate
initialize_dzhanibekov!(mech,
    angular_velocity=[15.0; 0.01; 0.0])
storage = simulate!(mech, 4.65,
    record=true,
    verbose=false)

# ## Visualizers
vis=visualizer()
render(vis)
visualize(mech, storage,
    vis=vis)
```

And voila! You should see something like this;

![dzhanibekov](./../../examples/animations/dzhanibekov.gif)
