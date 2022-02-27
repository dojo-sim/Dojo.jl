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
Ant
Atlas 
Box 
Cartpole
HalfCheetah
Hopper 
Pendulum
Quadruped 
RexHopper 
Walker 
```

## Simulation
```@docs
```

## Gradients
```@docs
```

## Visualization
```@docs
```


