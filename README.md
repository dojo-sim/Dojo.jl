# ConstrainedDynamics
[![Build Status](https://github.com/janbruedigam/ConstrainedDynamics.jl/workflows/CI/badge.svg)](https://github.com/janbruedigam/ConstrainedDynamics.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/janbruedigam/ConstrainedDynamics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/janbruedigam/ConstrainedDynamics.jl)
[![Dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://janbruedigam.github.io/ConstrainedDynamics.jl/dev/)

**ConstrainedDynamics** is a rigid body dynamics package written in Julia. In contrast to the RigidBodyDynamics package, it uses maximal coordinates instead of minimal coordinates to represent the states of a mechanism. This parametrization can be advantageous when simulating structures with additional explicit constraints. In terms of speed, ConstrainedDynamics and RigidBodyDynamics are roughly comparable, so for normal applications, both are valid options.

The package is largely built on StaticArrays and avoids allocations for improved performance. Convenience methods for setting up mechanical structures in maximal or minimal coordinates are provided and URDF parsing is also available. At the moment, there is a very basic documentation, but the provided examples should help with getting started with the package.

## Related Packages
* [RigidBodyDynamics](https://github.com/JuliaRobotics/RigidBodyDynamics.jl): Especially efficient for simulating unconstrained systems with single-degree-of-freedom joints
* [ConstrainedDynamicsVis](https://github.com/janbruedigam/ConstrainedDynamicsVis.jl): A package built on top of ConstrainedDynamics and MeshCat for visualization of simulations
* [ConstrainedControl](https://github.com/janbruedigam/ConstrainedControl.jl): Experimental implementation of Linear-Quadratic Regulation (LQR) in maximal coordinates

## Related Publications
[Linear-Time Variational Integrators in Maximal Coordinates](https://arxiv.org/abs/2002.11245) (accepted to WAFR 2020)

- [ ] we're not sure about this Bxmat, shouldn't this take into account some rot? trans?
```
    function g(mechanism, fric::Friction)
        body = getbody(mechanism, fric.parentid)
        x, v, q, ω = fullargssol(body.state)

        # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
        Bxmat = fric.Bx
```
