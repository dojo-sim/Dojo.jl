# Gradients

Convenient methods are provided for simulating an environment as well as computing the gradients of the dynamics with respect to the state and input (documentation for differentiation with respect to environment parameters coming soon!).

First, we will load an existing environment: 
```julia
env = get_environment(:cartpole, 
    representation=:maximal, 
    timestep=0.1,
    gravity=-9.81);
```

The environment can have `representation=:maximal` or `representation=:minimal`.

Next, we can inspect key environment dimensions:
```julia 
# dimensions
state_dim  = length(env.state)
jacobian_state_dim = size(env.dynamics_jacobian_state)
jacobian_input_dim = size(env.dynamics_jacobian_input)
```

The dynamics can be efficiently evaluted in place:
```julia
dynamics(next_state, env, current_state, input, parameters) 
```

and Jacobians can be similarly evaluated in place:
```julia
dynamics_jacobian_state(jacobian_state, env, current_state, input, parameters, 
    attitude_decompress=true)
dynamics_jacobian_input(jacobian_input, env, current_state, input, parameters, 
    attitude_decompress=true)
```
Note that the keyword argument `attitude_decompress` will determine whether the gradients account for special differentiation rules for quaternions (see [Planning with Attitude](https://roboticexplorationlab.org/papers/planning_with_attitude.pdf) for additional details). Simply, `attitude_decompress = true` will return Jacobians of standard dimensions without taking into account quaternion specialness.

For examples using this interface, see [Cart-pole](https://github.com/dojo-sim/Dojo.jl/blob/main/examples/trajectory_optimization/cartpole_max.jl).

