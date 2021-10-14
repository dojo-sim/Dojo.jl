# Mechanism from Code

Besides parsing a URDF file, you can also code up your mechanism.

```@example
## Dynamics simulation
using ConstrainedDynamics

# Origin
origin = Origin{Float64}()

# Links
link1 = Box(0.1, 0.1, 1, 1) # Box with length-x, length-y, length-z, and mass
link2 = Box(0.1, 0.1, 1, 1)

links = [link1;link2]

# Constraints
axis = [1;0;0] # rotation axis
connection1 = [0;0;0.5] # connection points in body frame
connection2 = [0;0;-0.5]
joint01 = EqualityConstraint(Revolute(origin, link1, axis; p2=connection1))
joint12 = EqualityConstraint(Revolute(link1, link2, axis; p1=connection2, p2=connection1))

constraints = [joint01;joint12]

# Mechanism
mech = Mechanism(origin, links, constraints)

# Set feasible initial configuration
setPosition!(origin, link1, p2 = connection1) # set link1 relative to origin
setPosition!(link1, link2, p1 = connection2, p2 = connection1, Δq = UnitQuaternion(RotX(pi/2))) # set link2 relative to link1 with an offset angle of pi/2

# Simulate
storage = simulate!(mech, 10, record = true) # simulate mech for 10 seconds and record the result (will be stored in storage)

## Visualization
using ConstrainedDynamicsVis 

visualize(mech, storage) # visualize mech with the data in storage.
 
```

