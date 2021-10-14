# Interface

## Initialization
The unused kwargs of initialization functions are set to zero.

```@autodocs
Modules = [ConstrainedDynamics]
Order   = [:function]
Pages   = ["initialize.jl"]
```

## Indexing
The following functions can be used to get bodies and constraints from a mechanism.

```@autodocs
Modules = [ConstrainedDynamics]
Order   = [:function]
Pages   = ["mechanism_functions.jl"]
```

## Minimal Coordinates
The state of a mechanism ist stored in maximal coordinates, i.e. the position and orientation of each body. Minimal coordinates (generalized, joint coordinates) can be set and retrived from constraints with the following functions.

```@autodocs
Modules = [ConstrainedDynamics]
Order   = [:function]
Pages   = ["equalityconstraint.jl"]
```
