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
maximal_dimension 
minimal_dimension
input_dimension
zero_velocity!
root_to_leaves_ordering
set_floating_base
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
Floating
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
contact_constraint
contact_location
get_sdf 
SphereHalfSpaceCollision
SphereSphereCollision
SphereCapsuleCollision
SphereBoxCollision
StringCollision
contact_normal 
contact_tangent
```

### Representations
```@docs
State
minimal_to_maximal
maximal_to_minimal
```

### Mechanics 
```@docs 
mechanical_energy 
kinetic_energy 
potential_energy 
momentum
```

### Graph
```@docs
System
Entry
full_matrix 
full_vector
```

## Environments
```@docs
Environment
Ant
Atlas
Block
Cartpole
HalfCheetah
Hopper
Pendulum
Quadruped
RaibertHopper
RexHopper
Walker
get_environment
step
get_observation
cost
is_done
reset
dynamics
dynamics_jacobian_state 
dynamics_jacobian_input
Space
BoxSpace
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
maximal_to_minimal_jacobian
minimal_to_maximal_jacobian
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
