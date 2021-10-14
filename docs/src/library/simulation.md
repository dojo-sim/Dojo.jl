# Simulation

## Simulation Methods

```@autodocs
Modules = [ConstrainedDynamics]
Order   = [:function]
Pages   = ["simulate.jl"]
```

## `Storage`

```@docs
Storage
```

## Control Function and `Controller`
Users can pass a control function to `simulate!`. This control function must be of the form
```julia
control_function(mechanism, timestep)
```
If such a control function is passed to `simulate!`, at each time step, `control_function(mechanism, timestep)` will be called. 
This way, the user has access to all attributes of the `mechanism` and the time step.

It is also possible to pass in a controller object:

```@docs
Controller
```

## Visualization
For details, see the `ConstrainedDynamicsVis` package.

Basic visualization can be performed with:
```julia
using ConstrainedDynamicsVis

visualize(mechanism, storage)
```