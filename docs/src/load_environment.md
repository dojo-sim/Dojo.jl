# Load Existing Environment

Dojo includes a number of existing environments: 

- [`Dojo.Ant`](@ref)
- [`Dojo.Atlas`](@ref)
- [`Dojo.Cartpole`](@ref)
- [`Dojo.HalfCheetah`](@ref)
- [`Dojo.Hopper`](@ref)
- [`Dojo.Pendulum`](@ref)
- [`Dojo.Quadruped`](@ref)
- [`Dojo.RaibertHopper`](@ref)
- [`Dojo.RexHopper`](@ref)
- [`Dojo.Walker`](@ref)

Specific environments can be instantiated, for example [`Dojo.Atlas`](@ref):

```julia 
env = get_environment(:atlas)
```