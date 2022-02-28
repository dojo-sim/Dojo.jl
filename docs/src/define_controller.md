# Defining a Controller
Here, we explain how to write a controller and simulate its effect on a dynamical system
i.e. a [`Mechanism`](@ref).
We focus on a simple pendulum swing-up.

Load Dojo and use the pendulum mechanism with desired simulation time step, desired gravity and desired damping at the joint.

```julia
using Dojo

mechanism = get_mechanism(:pendulum,
    timestep=0.01,
    gravity=-9.81,
    damper=5.0)
```

Define the controller. This is a method that takes 2 input arguments:
- a [`Mechanism`](@ref),
- an integer `k` indicating the current simulation step.
The controller computes the control inputs based on the current state `x`, the goal state `x_goal` and a proportional gain `K`.


```julia
function controller!(mechanism, k)
    ## Target state
    x_goal = [1.0 * π; 0.0]

    ## Current state
    x = get_minimal_state(mechanism)

    ## Gains
    K = [5.0 0.5] * 0.1

    # Control inputs
    u =  -K * (x - x_goal)
    set_input!(mechanism, u)
end
```

We initialize the pendulum in the lower position, and we will let the controller perform the swing-up movement.
```julia
initialize!(mechanism, :pendulum,
    angle=0.0 * π,
    angular_velocity=0.0);
```

We simulate the system for 2 seconds using the `controller!`.
```julia
storage = simulate!(mechanism, 2.0, controller!,
    record=true,
    verbose=true);
```

We visualize the results.
```julia
vis = Visualizer()
open(vis)
visualize(mechanism, storage, vis=vis);
```

You should get something like this,
```@raw html
<img src="./../../examples/animations/pendulum.gif" width="300"/>
```
