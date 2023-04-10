# Using an Existing Mechanism
The following code uses a function defined in `DojoEnvironments` to create a pendulum `Mechanism`. As before, the mechanism consists of an `Origin`, a `body`, and a `JointConstraint`. You can, of course, use these existing templates as a starting point for your own `Mechanism`s.

```julia
# ### Setup
using Dojo
using DojoEnvironments

# ### Get mechanism (check DojoEnvironment/mechanisms files for kwargs)
mechanism = get_mechanism(:pendulum; timestep=0.02, length=0.75) 
```
