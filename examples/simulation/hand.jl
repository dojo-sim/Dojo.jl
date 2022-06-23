using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo
using DojoEnvironments

# ## Mechanism
timestep = 0.01 
gravity = [0.0; 0.0; 0.0] 
spring = 0.0 
damper = 0.0
T = Float64

# Parameters
palm_mass = 1.0
palm_length = 1.0 
palm_radius = 0.5

segment_mass = 0.1
segment_length = 0.5
segment_radius = 0.2

rod_mass = 1.0
rod_length = 3.0 
rod_radius = 0.2

# Links
origin = Origin{T}()

palm = Capsule(palm_radius, palm_length, palm_mass,
    orientation_offset=RotX(0.5 * π),
    color=RGBA(0.7, 0.7, 0.7, 0.5))

# thumb
thumb_a = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

thumb_b = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

thumb_c = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

thumb = [thumb_a, thumb_b, thumb_c]

# index
index_a = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

index_b = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

index_c = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

index = [index_a, index_b, index_c]

# middle
middle_a = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

middle_b = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

middle_c = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

middle = [middle_a, middle_b, middle_c]

# ring
ring_a = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

ring_b = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

ring_c = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

ring = [ring_a, ring_b, ring_c]

# pinky
pinky_a = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

pinky_b = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

pinky_c = Capsule(segment_radius, segment_length, segment_mass, 
    orientation_offset=RotX(0.0 * π), 
    color=RGBA(0.7, 0.7, 0.7, 0.5))

pinky = [pinky_a, pinky_b, pinky_c]

rod = Capsule(rod_radius, rod_length, rod_mass, 
    orientation_offset=RotX(0.5 * π), 
    color=RGBA(1.0, 0.0, 0.0, 1.0))

links = [palm, thumb..., index..., middle..., ring..., pinky..., rod]

# Joints
origin_palm_joint = JointConstraint(Floating(origin, palm))

palm_rod_joint = JointConstraint(Floating(palm, rod))

# thumb
palm_thumb_a_joint = JointConstraint(Revolute(palm, thumb_a, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.5 * palm_length; -palm_radius], 
    child_vertex=[0.0; 0.0; 0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

thumb_a_thumb_b_joint = JointConstraint(Revolute(thumb_a, thumb_b, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; -0.5 * segment_length], 
    child_vertex=[0.0; 0.0; 0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

thumb_b_thumb_c_joint = JointConstraint(Revolute(thumb_b, thumb_c, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; -0.5 * segment_length], 
    child_vertex=[0.0; 0.0; 0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

thumb_joints = [palm_thumb_a_joint, thumb_a_thumb_b_joint, thumb_b_thumb_c_joint]

# index
palm_index_a_joint = JointConstraint(Revolute(palm, index_a, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.75 * palm_length; palm_radius], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

index_a_index_b_joint = JointConstraint(Revolute(index_a, index_b, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

index_b_index_c_joint = JointConstraint(Revolute(index_b, index_c, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

index_joints = [palm_index_a_joint, index_a_index_b_joint, index_b_index_c_joint]

# middle
palm_middle_a_joint = JointConstraint(Revolute(palm, middle_a, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.25 * palm_length; palm_radius], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

middle_a_middle_b_joint = JointConstraint(Revolute(middle_a, middle_b, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

middle_b_middle_c_joint = JointConstraint(Revolute(middle_b, middle_c, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

middle_joints = [palm_middle_a_joint, middle_a_middle_b_joint, middle_b_middle_c_joint]

# ring
palm_ring_a_joint = JointConstraint(Revolute(palm, ring_a, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; -0.25 * palm_length; palm_radius], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

ring_a_ring_b_joint = JointConstraint(Revolute(ring_a, ring_b, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

ring_b_ring_c_joint = JointConstraint(Revolute(ring_b, ring_c, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

ring_joints = [palm_ring_a_joint, ring_a_ring_b_joint, ring_b_ring_c_joint]

# pinky
palm_pinky_a_joint = JointConstraint(Revolute(palm, pinky_a, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; -0.75 * palm_length; palm_radius], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

pinky_a_pinky_b_joint = JointConstraint(Revolute(pinky_a, pinky_b, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

pinky_b_pinky_c_joint = JointConstraint(Revolute(pinky_b, pinky_c, [0.0; 1.0; 0.0]; 
    parent_vertex=[0.0; 0.0; 0.5 * segment_length], 
    child_vertex=[0.0; 0.0; -0.5 * segment_length],
    spring=spring, 
    damper=damper
    )
)

pinky_joints = [palm_pinky_a_joint, pinky_a_pinky_b_joint, pinky_b_pinky_c_joint]

joints = [origin_palm_joint, thumb_joints..., index_joints..., middle_joints..., ring_joints..., pinky_joints..., palm_rod_joint] 

# Mechanism
hand_mechanism = Mechanism(origin, links, joints;
    gravity, 
    timestep
)

# Initialize 

# thumb
# hand_mechanism.bodies[2].state.x2 = [0.0; 0.5 * palm_length; -1.0 * (0.5 * segment_length + palm_radius)]
# hand_mechanism.bodies[3].state.x2 = [0.0; 0.5 * palm_length; -1.0 * (0.5 * segment_length + palm_radius + 0.5 * segment_length)]
# hand_mechanism.bodies[4].state.x2 = [0.0; 0.5 * palm_length; -1.0 * (0.5 * segment_length + palm_radius + 1.0 * segment_length)]

# # index
# hand_mechanism.bodies[5].state.x2 = [0.0; 0.75 * palm_length; 1.0 * (0.5 * segment_length + palm_radius)]
# hand_mechanism.bodies[6].state.x2 = [0.0; 0.75 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 0.5 * segment_length)]
# hand_mechanism.bodies[7].state.x2 = [0.0; 0.75 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 1.0 * segment_length)]

# # middle
# hand_mechanism.bodies[8].state.x2 = [0.0; 0.25 * palm_length; 1.0 * (0.5 * segment_length + palm_radius)]
# hand_mechanism.bodies[9].state.x2 = [0.0; 0.25 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 0.5 * segment_length)]
# hand_mechanism.bodies[10].state.x2 = [0.0; 0.25 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 1.0 * segment_length)]

# # middle
# hand_mechanism.bodies[11].state.x2 = [0.0; -0.25 * palm_length; 1.0 * (0.5 * segment_length + palm_radius)]
# hand_mechanism.bodies[12].state.x2 = [0.0; -0.25 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 0.5 * segment_length)]
# hand_mechanism.bodies[13].state.x2 = [0.0; -0.25 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 1.0 * segment_length)]

# # pinky 
# hand_mechanism.bodies[14].state.x2 = [0.0; -0.75 * palm_length; 1.0 * (0.5 * segment_length + palm_radius)]
# hand_mechanism.bodies[15].state.x2 = [0.0; -0.75 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 0.5 * segment_length)]
# hand_mechanism.bodies[16].state.x2 = [0.0; -0.75 * palm_length; 1.0 * (0.5 * segment_length + palm_radius + 1.0 * segment_length)]

# rod 
hand_mechanism.bodies[17].state.x2 = [palm_radius + rod_radius + 0.0; 0.0; 0.0]

state_minimal = get_minimal_state(hand_mechanism)

state_minimal[25] = 0.5 * π
state_minimal[27] = 0.1 * π
state_minimal[29] = 0.25 * π

state_minimal[31] = 0.5 * π
state_minimal[33] = 0.1 * π
state_minimal[35] = 0.25 * π

state_minimal[37] = 0.5 * π
state_minimal[39] = 0.1 * π
state_minimal[41] = 0.25 * π

state_minimal[43] = 0.5 * π
state_minimal[45] = 0.1 * π
state_minimal[47] = 0.25 * π

state_minimal[49] = -0.3 * π
state_minimal[51] = -0.1 * π
state_minimal[53] = -0.25 * π

set_minimal_state!(hand_mechanism, state_minimal)

# Simulate
storage = simulate!(hand_mechanism, timestep, 
    record=true, 
    opts=SolverOptions(
        rtol=1.0e-4, 
        btol=1.0e-4, 
        verbose=true
    )
)

# Visualize
# vis = Visualizer()
# render(vis)
visualize(hand_mechanism, storage, 
    vis=vis, 
)

setvisible!(vis[:floor], false)
