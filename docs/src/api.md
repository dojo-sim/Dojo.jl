# API Documentation

Docstrings for Dojo.jl interface members can be [accessed through Julia's built-in documentation system](https://docs.julialang.org/en/v1/manual/documentation/index.html#Accessing-Documentation-1) or in the list below.

```@meta
CurrentModule = Dojo
```

## Contents

```@contents
Pages = ["api.md"]
```

## Index

```@index
Pages = ["api.md"]
```

## Mechanism

```@docs
Mechanism
get_mechanism 
initialize!
get_node 
get_body
get_joint 
get_contact 
get_maximal_state 
get_next_state
get_minimal_state
set_maximal_state!
set_minimal_state!
set_input! 
input_dimension
```

### Nodes 
```@docs
Node
Body
Origin
Constraint
Shape
EmptyShape 
Mesh
Box
Cylinder 
Capsule
Sphere
Pyramid
Shapes
```

### Joints
```@docs
Joint 
Rotational 
Translational
JointConstraint 
Fixed
Prismatic
Planar
FixedOrientation
Revolute
Cylindrical
PlanarAxis
FreeRevolute
Orbital
PrismaticOrbital
PlanarOrbital
FreeOrbital
Spherical
CylindricalFree
PlanarFree
```

### Contacts
```@docs
Contact 
ImpactContact
LinearContact 
NonlinearContact
ContactConstraint
```

### Representations
```@docs
State
```

### Graph 
```@docs 
System
Entry
```

## Environments 
```@docs
Environment
Space 
BoxSpace
Ant
Atlas 
Block 
Cartpole
HalfCheetah
Hopper 
Pendulum
Quadruped 
RexHopper 
Walker 
get_environment 
step
get_observation 
cost 
is_done
reset
```

## Simulate
```@docs
    Storage 
    step!
    simulate! 
```

## Gradients
```@docs
get_maximal_gradients! 
get_minimal_gradients!
```

## Solver 
```@docs
SolverOptions
mehrotra!
```

## Visualization
```@docs
visualize
build_robot 
set_robot
set_camera! 
set_light! 
set_surface! 
set_floor!
```


