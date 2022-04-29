# Simulation

The results of a simulation step can be differentiated with respect to current state and control input. (Documentation for differentiation with respect to mechcanism parameters coming soon!).

See [environment simulation with gradients](environment_simulation_gradients.md) for more convenient methods for simulating an [`Dojo.Environment`](@ref) and computing gradients.

```julia
using Dojo

# load mechanism
mechanism = get_mechanism(:cartpole,
    timestep=0.1,
    gravity=-9.81)

# get current maximal state
maximal_state = get_maximal_state(mechanism)

# random input
input = randn(input_dimension(mechanism))

# maximal state and input jacobians of the mechanism dynamics
jacobian_maximal_state, jacobian_input = get_maximal_gradients!(mechanism, maximal_state, input)

# minimal state and input jacobians of the mechanism dynamics
jacobian_minimal_state, jacobian_input = get_minimal_gradients!(mechanism, maximal_state, input)
```

