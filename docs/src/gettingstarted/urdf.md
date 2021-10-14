# Mechanism from URDF

The easiest way to run a simulation is to parse a URDF file.

```@example
## Dynamics simulation
using ConstrainedDynamics

path = "path_to_file/doublependulum.urdf"
mech = Mechanism(path; Δt = 0.01)
storage = simulate!(mech, 10, record = true) # simulate mech for 10 seconds and record the result (will be stored in storage)

## Visualization
using ConstrainedDynamicsVis 

visualize(mech, storage) # visualize mech with the data in storage.
 
```

