# Using an Existing Environment
The following code uses a function defined in `DojoEnvironments` to create a pendulum `Environment`. This `Environment` is a wrapper around the [`Mechanism`](@ref) for easy interfacing with other packages. As with [`Mechanism`](@ref)s, you can use the existing templates as a starting point for your own `Environment`s.

```julia
# ### Setup
using Dojo
using DojoEnvironments

# ### Get environment (check DojoEnvironment/environments files for kwargs)
environment = get_environment(:pendulum; timestep=0.01, horizon=200) 
```
