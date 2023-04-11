# Directly Build a Mechanism
The following code builds a pendulum [`Mechanism`](@ref), consisting of an [`Origin`](@ref), a [`Body`](@ref), and a [`JointConstraint`](@ref). Most [`Mechanism`](@ref)s consist of these three components, sometimes supplemented by [`ContactConstraint`](@ref)s. 

```julia
# ### Setup
using Dojo

# ### Parameters
radius = 0.1
length = 1
mass = 1
rotation_axis = [1;0;0] 
connection = [0;0;length/2]

# ### Mechanism components
origin = Origin()
body = Cylinder(radius, length, mass)
joint = JointConstraint(Revolute(origin, body, rotation_axis; child_vertex=connection))

# ### Construct Mechanism
mechanism = Mechanism(origin, [body], [joint])
```
