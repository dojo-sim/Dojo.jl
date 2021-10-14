# ConstrainedDynamics

**ConstrainedDynamics** is a rigid body dynamics package written in Julia. It uses maximal coordinates instead of minimal coordinates to represent the states of a mechanism. This parametrization can be advantageous when simulating structures with additional explicit constraints.

The package is largely built on StaticArrays and avoids allocations for improved performance. Convenience methods for setting up mechanical structures in maximal or minimal coordinates are provided and URDF parsing is also available. The provided examples should help with getting started with the package. 


## Contents

```@contents
Pages = [
  "body.md",
  "constraint.md",
  "mechanism.md",
  "interface.md",
  "simulation.md"]
Depth = 2
```