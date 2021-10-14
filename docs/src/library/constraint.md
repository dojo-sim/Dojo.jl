# Constraints

## `EqualityConstraint`

```@docs
EqualityConstraint
```

## Joint Prototypes
The functions below can be used to create joints between two bodies. Note that body1 can also be the origin. The following arguments are available for most joints:
* `axis`:     The rotation axis or plane axis.
* `p1`:       The joint connection point for body1 in body1's frame
* `p2`:       The joint connection point for body2 in body2's frame
* `qoffsett`: The orientation offset of body2 relative to body1

```@autodocs
Modules = [ConstrainedDynamics]
Order   = [:function]
Pages   = ["prototypes.jl"]
```
