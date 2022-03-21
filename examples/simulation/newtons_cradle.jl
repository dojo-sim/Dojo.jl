timestep=0.1
gravity=-9.81

T=Float64

# Parameters
height = 1.0
num_rods = 3
rod_separation = 0.1

mass_connector = 1.0 
length_connector = (num_rods - 1) * rod_separation
radius_connector = 0.025

mass_rod = 0.1 
length_rod = 0.5
radius_rod = 0.01

# Bodies
origin = Origin{T}()
connector = Capsule(radius_connector, length_connector, mass_connector, 
    axis_offset=Quaternion(0.707107, 0.0, 0.707107, 0.0),
    name=:connector)

rods = [Capsule(radius_rod, length_rod, mass_rod, name=Symbol("rod_$i")) for i = 1:num_rods]

# Constraints
joint_axis = [0.0; 1.0; 0.0]

joint_world_connector = JointConstraint(Fixed(origin, connector, parent_vertex=[0.0; 0.0; height]))

joint_connector_rods = [JointConstraint(Revolute(connector, rod, joint_axis; parent_vertex=[(i - 1) * rod_separation - 0.5 * rod_separation * (num_rods - 1); 0.0; 0.0], child_vertex=[0.0; 0.0; 0.5 * length_rod]), name=Symbol("connector_rod_$i")) for (i, rod) in enumerate(rods)]

bodies = [connector, rods...]
joints = [joint_world_connector, joint_connector_rods...]

mech = Mechanism(origin, bodies, joints,
    gravity=gravity,
    timestep=timestep)

set_minimal_state!(mech, [-0.1; 0.0; 0.0; 0.0; 0.0; 0.0])

storage = simulate!(mech, 1.0, 
    record=true, 
    opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6, verbose=true))

visualize(mech, storage, 
    vis=vis, show_contact=true)

# ## Visualize
# vis = Visualizer()
# render(vis)
# build_robot(mech, vis=vis)
# set_robot(vis, mech, get_maximal_state(mech))



