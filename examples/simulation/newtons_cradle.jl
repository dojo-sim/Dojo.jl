timestep = 1.0e-3
gravity = -9.81

T = Float64

# Parameters
height = 1.0
num_balls = 2
rod_separation = 0.1

mass_connector = 1.0
length_connector = (num_balls - 1) * rod_separation
radius_connector = 0.025

mass_rod = 0.001
length_rod = 0.5
radius_rod = 0.01

mass_ball = 1.0
radius_ball = 0.5 * rod_separation

# Bodies
origin = Origin{T}()
<<<<<<< HEAD
connector = Capsule(radius_connector, length_connector, mass_connector,
    axis_offset=Quaternion(0.707107, 0.707107, 0.0,  0.0),
=======
connector = Capsule(radius_connector, length_connector, mass_connector, 
    orientation_offset=Quaternion(0.707107, 0.707107, 0.0,  0.0),
>>>>>>> 33ccdde8d7f74c4ea3c3bffdebcc6ed65959a5be
    name=:connector)

rods = [Capsule(radius_rod, length_rod, mass_rod,
    name=Symbol("rod_$i")) for i = 1:num_balls]
balls = [Sphere(radius_ball, mass_ball,
    name=Symbol("ball_$i")) for i = 1:num_balls]

# Joints
joint_axis = [1.0; 0.0; 0.0]
damper = 0.0
joint_world_connector = JointConstraint(Fixed(origin, connector,
    parent_vertex=[0.0; 0.0; height]))

joint_connector_rods = [JointConstraint(Revolute(connector, rod, joint_axis;
    damper=damper,
    parent_vertex=[0.0; (i - 1) * rod_separation - 0.5 * rod_separation * (num_balls - 1);  0.0],
    child_vertex=[0.0; 0.0; 0.5 * length_rod]),
    name=Symbol("connector_rod_$i")) for (i, rod) in enumerate(rods)]

# i = 1
# fixed_joint = JointConstraint(Fixed(connector, rods[i];
#         parent_vertex=[0.0; (i - 1) * rod_separation - 0.5 * rod_separation * (num_balls - 1);  0.0],
#         child_vertex=[0.0; 0.0; 0.5 * length_rod]),
#         name=Symbol("connector_rod_$i"))

# joint_connector_rods = [fixed_joint, joint_connector_rods[2]]


joint_rods_balls = [JointConstraint(Fixed(rods[i], balls[i],
    parent_vertex=[0.0; 0.0; -0.5 * length_rod],
    child_vertex=[0.0; 0.0; 0.0])) for i = 1:num_balls]

# Contacts
friction_coefficient = 0.5

collision0 = SphereSphereCollision{Float64,2,3,6}(
    [0.0; 0.0; 0.0],
    [0.0; 0.0; 0.0],
    radius_ball,
    radius_ball)

friction_parameterization = SA{Float64}[
        1.0  0.0
        0.0  1.0
]

body_body = NonlinearContact{Float64,8}(friction_coefficient, friction_parameterization, collision0)

# collision0 = SphereSphereCollision{Float64,0,3,0}(
#     [0.0; 0.0; 0.0],
#     [0.0; 0.0; 0.0],
#     radius_ball,
#     radius_ball)

# friction_parameterization = szeros(Float64, 0, 2)

# body_body = ImpactContact{Float64,8}(friction_parameterization, collision0)

bodies = [connector, rods..., balls...]
joints = [joint_world_connector, joint_connector_rods..., joint_rods_balls...]
contacts = [RigidContactConstraint((body_body, bodies[i + 1 + num_balls].id, bodies[i + 1 + num_balls + 1].id),
    name=Symbol("body_body_$i")) for i = 1:(num_balls-1)]

mech = Mechanism(origin, bodies, joints, contacts,
    gravity=gravity,
    timestep=timestep)

minimal_dimension(mech)
set_minimal_state!(mech, [0.25; 0.0; 0.0; 0.0; 0.0; 0.0])

@show adjacency_matrix(mech.joints, mech.bodies, mech.contacts)

distance(mech.contacts[1].model.collision,
    mech.bodies[1 + num_balls + 1].state.x2, mech.bodies[1 + num_balls + 1].state.q2,
    mech.bodies[1 + num_balls + 2].state.x2, mech.bodies[1 + num_balls + 2].state.q2,)

contact_point(:parent, mech.contacts[1].model.collision,
    mech.bodies[1 + num_balls + 1].state.x2, mech.bodies[1 + num_balls + 1].state.q2,
    mech.bodies[1 + num_balls + 2].state.x2, mech.bodies[1 + num_balls + 2].state.q2,)

contact_point(:child, mech.contacts[1].model.collision,
    mech.bodies[1 + num_balls + 1].state.x2, mech.bodies[1 + num_balls + 1].state.q2,
    mech.bodies[1 + num_balls + 2].state.x2, mech.bodies[1 + num_balls + 2].state.q2,)

contact_normal(mech.contacts[1].model.collision,
    mech.bodies[1 + num_balls + 1].state.x2, mech.bodies[1 + num_balls + 1].state.q2,
    mech.bodies[1 + num_balls + 2].state.x2, mech.bodies[1 + num_balls + 2].state.q2,)

contact_tangent(mech.contacts[1].model.collision,
    mech.bodies[1 + num_balls + 1].state.x2, mech.bodies[1 + num_balls + 1].state.q2,
    mech.bodies[1 + num_balls + 2].state.x2, mech.bodies[1 + num_balls + 2].state.q2,)

storage = simulate!(mech, 1.0,
    record=true,
    opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6, verbose=true))

vis = Visualizer()
render(vis)
visualize(mech, storage,
    vis=vis)
open(vis)

# ## Visualize
# vis = Visualizer()
# render(vis)
build_robot(mech, vis=vis)
set_robot(vis, mech, get_maximal_state(mech), show_contact=false)

convert_frames_to_video_and_gif("cradle")
