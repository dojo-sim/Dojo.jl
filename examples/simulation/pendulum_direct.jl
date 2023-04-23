# ### Setup
# PKG_SETUP
using Dojo

# ### Parameters
radius = 0.1
link_length = 1
mass = 1
rotation_axis = [1;0;0] 
connection = [0;0;link_length/2]

# ### Mechanism components
origin = Origin()
body = Cylinder(radius, link_length, mass)
joint = JointConstraint(Revolute(origin, body, rotation_axis; child_vertex=connection))

# ### Construct Mechanism
mechanism = Mechanism(origin, [body], [joint])

# ### Set state
set_minimal_coordinates!(mechanism, joint, [pi/4])
set_minimal_velocities!(mechanism, joint, [0.2])

# ### Simulate
storage = simulate!(mechanism, 5.0, record=true)

# ### Visualize
vis = visualize(mechanism, storage; visualize_floor=false)
render(vis)